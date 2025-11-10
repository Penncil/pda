# Copyright 2021 Penn Computing Inference Learning (PennCIL) lab
#       https://penncil.med.upenn.edu/team/
# This file is part of pda
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# https://style.tidyverse.org/functions.html#naming
# https://ohdsi.github.io/Hades/codeStyle.html#OHDSI_code_style_for_R

COLA.steps <- c('initialize','estimate')
COLA.family <- 'binomial'

# One-shot lossless generalized linear model
# 3 models: 
# - COLA-GLM
# - COLA-GLM-H
# - COLA-GLMM

# Helper 
`%||%` <- function(a, b) if (is.null(a)) b else a


#' Construct binary covariate pattern matrix
#'
#' Builds the full design grid of binary fixed-effect patterns used by COLA
#' to aggregate counts and outcome sums. When `intercept = TRUE`, an `Intercept`
#' column of ones is included and all other variables are expanded over \{0,1\}.
#'
#' @param x_names Character vector of fixed-effect names. If `intercept = TRUE`,
#'   it may include `"Intercept"`; otherwise it must not.
#' @param intercept Logical; include a fixed intercept column. Default: `TRUE`.
#'
#' @return A tibble of all binary patterns over `x_names` (with/without `Intercept`),
#'   one row per pattern. 
#' @keywords internal
make_patterns <- function(x_names, intercept = TRUE) {
  x_names <- as.character(x_names %||% character())
  
  if (intercept) {
    bin_names <- setdiff(x_names, "Intercept")
    Xbin <- if (length(bin_names)) {
      expand.grid(rep(list(c(0L, 1L)), length(bin_names)))
    } else {
      data.frame(dummy = integer())
    }
    if (ncol(Xbin) == 0) Xbin <- data.frame()
    names(Xbin) <- bin_names
    tibble::as_tibble(cbind(Intercept = 1L, Xbin))
  } else {
    if ("Intercept" %in% x_names) stop("When intercept = FALSE, do not include 'Intercept' in x_names.")
    if (length(x_names) == 0) return(tibble::as_tibble(data.frame()))
    Xbin <- expand.grid(rep(list(c(0L, 1L)), length(x_names)))
    names(Xbin) <- x_names
    tibble::as_tibble(Xbin)
  }
}



#' One-shot site summaries for COLA-GLMM
#'
#' Produces the **lossless**, pattern-level sufficient statistics for a single site:
#' pattern counts `Ck`, outcome sums `Sk = Σy`, squared sums `S2k = Σy^2`, and the
#' corresponding pattern matrix `X0`. Works for both binomial and Poisson outcomes
#' (for Bernoulli, `S2k == Sk`).
#'
#' @param df_site Data frame for one site. Must include outcome column named `y`
#'   and the fixed-effect covariates in `x_names`. If `intercept = TRUE`, the
#'   function will add an `Intercept` column when missing.
#' @param x_names Character vector of fixed-effect names (binary covariates; may
#'   include `"Intercept"` if `intercept = TRUE`).
#' @param intercept Logical; include a fixed intercept in the pattern matrix.
#'
#' @return A list with elements:
#' \itemize{
#'   \item `Ck`   (integer vector) pattern counts
#'   \item `Sk`   (numeric vector) sums of y per pattern
#'   \item `S2k`  (numeric vector) sums of y^2 per pattern
#'   \item `X0`   (matrix) pattern design matrix aligned to `Ck/Sk/S2k`
#' }
#' @examples
#' # df_site$y must exist; x_names are binary
#' # out <- generate_CSU_site(df_site, c("Intercept","age","sex"), intercept = TRUE)
#' @keywords internal
generate_CSU_site <- function(df_site,
                              x_names,
                              intercept = TRUE) {
  x_names <- as.character(x_names %||% character())
  
  if (intercept && !"Intercept" %in% x_names) {
    x_names <- c("Intercept", x_names)
  }
  
  X0 <- make_patterns(x_names, intercept = intercept)
  df_work <- df_site
  if (intercept && !"Intercept" %in% names(df_work)) {
    df_work <- dplyr::mutate(df_work, Intercept = 1L)
  }
  
  obs <- df_work |>
    dplyr::count(dplyr::across(dplyr::all_of(x_names)), name = "n") |>
    dplyr::left_join(
      df_work |>
        dplyr::group_by(dplyr::across(dplyr::all_of(x_names))) |>
        dplyr::summarise(SY = sum(y), S2 = sum(y^2), .groups = "drop"),
      by = x_names
    )
  
  tb <- X0 |>
    dplyr::left_join(obs, by = x_names) |>
    dplyr::mutate(
      n  = tidyr::replace_na(n,  0L),
      SY = tidyr::replace_na(SY, 0L),
      S2 = tidyr::replace_na(S2, 0L)
    )
  
  list(
    Ck  = tb$n,
    Sk  = tb$SY,
    S2k = tb$S2,
    X0  = as.matrix(tb[, x_names, drop = FALSE])
  )
}


#' @useDynLib pda
#' @title COLA initialize
#' 
#' @usage COLA.initialize(ipdata,control,config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#' 
#' @references Qiong Wu, et al. (2025) COLA-GLM: Collaborative One-shot and Lossless Algorithms of Generalized Linear Models for Decentralized Observational Healthcare Data. 
#'    npj Digital Medicine. \cr
#'    Bingyu Zhang, et al (2025) A Lossless One-shot Distributed Algorithm for Addressing Heterogeneity in Multi-Site Generalized Linear Models.
#'    Journal of the American Medical Informatics Association (under revision). \cr
#'    Jiayi Tong, et al. (2025) Unlocking Efficiency in Real-world Collaborative Studies: A Multi-site International Study with Collaborative One-shot Lossless Algorithm for Generalized Linear Mixed Model.
#'    npj Digital Medicine. \cr
#'    
#' @return init
#' @keywords internal
COLA.initialize <- function(ipdata, control, config) {
  # --- COMMON ---
  ycol <- control$outcome %||% "outcome"
  vns  <- control$variables %||% setdiff(colnames(ipdata), ycol)
  intercept <- control$intercept %||% TRUE
  mixed <- control$mixed_effects %||% FALSE
  
  # --- UPDATED GLMM BRANCH ---
  if (isTRUE(mixed)) {
    df_site <- as.data.frame(ipdata[, c(vns, ycol), with = FALSE])
    names(df_site)[names(df_site) == ycol] <- "y"
    csu <- generate_CSU_site(df_site, x_names = vns, intercept = intercept) # get summary state step
    csu$x_names <- colnames(csu$X0)
    return(csu)
  }
  
  # --- GLM / GLM-H part ---
  Xmat <- ipdata[,!colnames(ipdata) %in% "outcome", with = FALSE]
  Y <- ipdata[,colnames(ipdata) %in% "outcome", with = FALSE]
  Xmat.tbl <- data.frame(Xmat)
  category_combinations <- expand.grid(lapply(Xmat.tbl, unique),
                                       stringsAsFactors = FALSE) %>% dplyr::arrange_all()
  colnames(Xmat.tbl) <- colnames(category_combinations)
  
  if (control$link == "canonical") {
    Xmat.tbl <- tibble::as_tibble(Xmat.tbl)
    cols <- colnames(Xmat.tbl)
    Xtable_initial <- Xmat.tbl %>%
      dplyr::group_by(across(everything())) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop")
    Xtable <- category_combinations %>%
      dplyr::left_join(Xtable_initial, by = cols) %>%
      as.data.frame()
    Xtable$n[which(is.na(Xtable$n))] = 0
    colnames(Xtable) <- c(colnames(Xmat),'n')
    if(is.numeric(control$cutoff)==TRUE){
      Xtable$n[Xtable$n>0 & Xtable$n <control$cutoff] <- rep(ceiling(control$cutoff/2))
    }
    
    if(control$mixed_effects==FALSE){ # glm
      SXY <- t(Y)%*%as.matrix(Xmat)
      init <- list(SXY = SXY, Xtable = Xtable)
    }
  }
  return(init)
}



#' COLA-GLMM 
#'
#' Fits a generalized linear mixed model with site-level random intercepts using
#' only one-shot per-site summaries \code{(Ck, Sk, S2k, X0)}. Each of the
#' iterations constructs weighted LMM summary statistics which are then solved
#' by \code{lmm.fit} (from DLMM), yielding updated fixed effects and random
#' intercepts until convergence.
#'
#' @param summary_by_site Named list of site summaries. Each element must contain
#'   \code{Ck}, \code{Sk}, \code{S2k}, and \code{X0} as returned by
#'   \code{generate_CSU_site}. The list names should be site IDs.
#' @param family Character; one of \code{"poisson"} or \code{"binomial"}
#'   (canonical links).
#' @param intercept Logical; whether the fixed-effect design includes an
#'   intercept (affects how \code{X0} was constructed). Default \code{TRUE}.
#' @param beta_init Optional named numeric vector of initial fixed effects.
#'   Defaults to zeros.
#' @param u_init Optional named numeric vector of initial site random effects
#'   (one per site). Defaults to zeros.
#' @param max_iter Integer maximum number of IRLS iterations. Default \code{50}.
#' @param tol Convergence tolerance on relative squared parameter change.
#'   Default \code{1e-6}.
#' @param verbose Logical; print iteration progress. Default \code{TRUE}.
#' @return A list with elements:
#'   \itemize{
#'     \item \code{beta}: named fixed-effect estimates
#'     \item \code{u}: named site random-intercept BLUPs
#'     \item \code{V}: variance component matrix for the random intercept
#'     \item \code{s2}: residual scale from the working LMM
#'     \item \code{iter}: number of iterations performed
#'     \item \code{SiXYZ_last}: last iteration's sufficient statistics (by site)
#'   }
#' @details
#' Uses canonical links: log for Poisson and logit for binomial. The fixed-effect
#' covariates in \code{X0} are assumed binary (plus optional \code{Intercept}).
#' For numerically extreme logits, a small weight floor is used internally.
#' Requires \code{lmm.fit} from \emph{dlmm.R} to be on the search path.
#' @examples
#' # fit <- cola_glmm(summary_by_site, family = "poisson")
#' @keywords internal
cola_glmm <- function(summary_by_site,
                      family    = "poisson",
                      intercept = TRUE,
                      beta_init = NULL, u_init = NULL,
                      max_iter  = 50, tol = 1e-6, verbose = TRUE) {
  
  sites <- names(summary_by_site)
  K     <- length(sites)
  
  X0 <- summary_by_site[[1]]$X0
  p  <- ncol(X0)
  x_names <- colnames(X0) %||% paste0("x", seq_len(p))
  
  beta <- if (is.null(beta_init)) setNames(rep(0, p), x_names) else beta_init
  u    <- if (is.null(u_init))    setNames(rep(0, K), sites)   else u_init
  
  # canonical links
  link_mu <- switch(family,
                    poisson  = function(eta) exp(eta),
                    binomial = function(eta) plogis(eta),
                    stop("family must be 'poisson' or 'binomial'"))
  
  w_fun <- switch(family,
                  poisson  = function(mu) mu,
                  binomial = function(mu) mu * (1 - mu))
  
  # site-level summaries to LMM-style blocks
  site_to_SiXYZ <- function(Ck, Sk, S2k, X0, u_k, beta) {
    n   <- as.numeric(Ck)
    sy  <- as.numeric(Sk)
    sy2 <- as.numeric(S2k)
    eta <- as.vector(X0 %*% beta) + u_k
    mu  <- link_mu(eta)
    w   <- pmax(w_fun(mu), 1e-12)
    
    nw     <- n * w
    sum_w  <- sum(nw)
    sum_wx <- colSums(X0 * nw)
    wz     <- nw * eta + (sy - n * mu)
    wz2    <- n * w * eta^2 + 2 * eta * (sy - n * mu) +
      (sy2 - 2 * mu * sy + n * mu^2) / w
    
    WX   <- X0 * as.numeric(nw)
    SiX  <- crossprod(WX, X0)
    SiXZ <- matrix(sum_wx, nrow = ncol(X0))
    SiZ  <- matrix(sum_w, nrow = 1, ncol = 1)
    SiXY <- matrix(colSums(X0 * as.numeric(wz)), nrow = ncol(X0))
    SiZY <- matrix(sum(wz), nrow = 1, ncol = 1)
    SiY  <- sum(wz2)
    
    list(SiX = SiX, SiXZ = SiXZ, SiXY = SiXY,
         SiZ = SiZ, SiZY = SiZY, SiY = SiY, ni = sum(n))
  }
  
  bu_prev <- c(beta, u)
  fit <- NULL
  
  for (it in seq_len(max_iter)) {
    if (verbose) message(sprintf("COLA-GLMM (%s) iteration %d", family, it))
    
    SiXYZ <- setNames(vector("list", K), sites)
    for (k in seq_len(K)) {
      ss <- summary_by_site[[k]]
      SiXYZ[[k]] <- site_to_SiXYZ(Ck = ss$Ck, Sk = ss$Sk, S2k = ss$S2k,
                                  X0 = ss$X0, u_k = u[[k]], beta = beta)
    }
    
    fit <- lmm.fit(
      pooled    = FALSE,
      reml      = TRUE,
      common.s2 = TRUE,
      SiXYZ     = SiXYZ,
      corstr    = "independence",
      verbose   = FALSE
    )
    
    beta <- setNames(as.numeric(fit$b), x_names)
    u    <- setNames(sapply(fit$ui, function(v) v[1]), sites)
    
    bu_new <- c(beta, u)
    relchg <- sum((bu_new - bu_prev)^2) / max(1e-12, sum(bu_prev^2))
    if (verbose) message(sprintf("  relative change = %.3e", relchg))
    bu_prev <- bu_new
    if (relchg < tol) break
  }
  
  list(beta = beta, u = u, V = fit$V, s2 = fit$s2, iter = it, SiXYZ_last = SiXYZ)
}



#' @useDynLib pda
#' @title PDA COLA estimation
#' 
#' @usage COLA.estimate(ipdata=NULL,control,config)
#' @param ipdata no need
#' @param control PDA control
#' @param config cloud configuration
#' 
#' @details COLA estimation: 
#' (1) COLA-GLM
#' (2) COLA-GLM-H
#' (3) COLA-GLMM
#' @return  list(est, se) 
#' @keywords internal
COLA.estimate <- function(ipdata=NULL, control, config) {
  K <- length(control$sites)
  
  if(control$mixed_effects==TRUE){
    # --- UPDATED GLMM BRANCH ---
    sites <- control$sites
    summary_by_site <- setNames(vector("list", length(sites)), sites)
    for (i in seq_along(sites)) {
      AD <- pdaGet(paste0(sites[i], "_initialize"), config)
      summary_by_site[[i]] <- AD
    }
    
    fam <- control$family %||% "poisson"
    intercept <- control$intercept %||% TRUE
    
    fit <- cola_glmm(summary_by_site = summary_by_site,
                     family    = fam,
                     intercept = intercept,
                     beta_init = control$beta_init %||% NULL,
                     u_init    = control$u_init %||% NULL,
                     max_iter  = control$max_iter %||% 50,
                     tol       = control$tol %||% 1e-6,
                     verbose   = control$verbose %||% TRUE)
    
    res <- list(est = fit$beta, u = fit$u, V = fit$V, s2 = fit$s2, iter = fit$iter)
    
  } else if(control$mixed_effects==FALSE){
    # --- GLM / GLM-H BRANCH ---
    if(control$heterogeneity==FALSE){
      # COLA-GLM
      AD <- pdaGet(paste0(control$sites[1],'_initialize'),config)
      SXY <- AD$SXY
      Xtable <- as.data.frame(matrix(unlist(AD$Xtable), ncol = (length(control$variables) + 2)))
      colnames(Xtable) <- c("intercept", control$variables, "n")
      Xcat <- as.matrix(Xtable[,colnames(Xtable)!='n'])
      counts <- Xtable$n
      for(site_i in control$sites[-1]){
        KSiteAD <- pdaGet(paste0(site_i,'_initialize'),config)
        SXY <- SXY+KSiteAD$SXY 
        Xtable <- as.data.frame(matrix(unlist(KSiteAD$Xtable), ncol = (length(control$variables) + 2)))
        colnames(Xtable) <- c("intercept", control$variables, "n")
        counts <- counts + Xtable$n
      }
    } else {
      # COLA-GLM-H
      AD <- pdaGet(paste0(control$sites[1],'_initialize'),config)
      SXY.intercept <- AD$SXY[1]
      SXY.cov <- AD$SXY[-1]
      Xtable <-as.data.frame(matrix(unlist(AD$Xtable), ncol = (length(control$variables) + 2)))
      colnames(Xtable) <- c("intercept", control$variables, "n") 
      Xcat0 <- Xtable[,colnames(Xtable)!='n']
      Xcat <- Xcat0[,-1]
      counts <- Xtable$n
      
      for(site_i in control$sites[-1]){
        KSiteAD <- pdaGet(paste0(site_i,'_initialize'),config)
        SXY.intercept <- c(SXY.intercept,KSiteAD$SXY[1])
        SXY.cov <- SXY.cov+KSiteAD$SXY[-1] 
        Xtable <- as.data.frame(matrix(unlist(KSiteAD$Xtable), ncol = (length(control$variables) + 2)))
        colnames(Xtable) <- c("intercept", control$variables, "n")
        Xcat0 <- Xtable[,colnames(Xtable)!='n']
        Xcat_tmp <- Xcat0[,-1]
        Xcat <- rbind(Xcat,Xcat_tmp)
        counts <- c(counts, Xtable$n)
      }
      
      SXY <- c(SXY.intercept,SXY.cov)
      SiteID <- rep(1:K, each = nrow(Xcat0))
      new.siteID <- sapply(1:K,function(i) ifelse(SiteID==i,1,0))
      colnames(new.siteID) <- paste0("Site", 1:K)
      Xcat <- cbind(new.siteID,Xcat)
      Xcat <- as.matrix(Xcat)
    }
    
    logLik_AD <- function(beta){
      if(control$family == 'binomial'){
        -(sum(SXY*beta)-sum(log1p(exp(Xcat%*%c(beta)))*counts))/sum(counts)
      } else if (control$family == 'poisson'){
        lambda <- exp(Xcat %*% beta)
        -(sum(SXY * beta) - sum(lambda * counts)) / sum(counts)
      }
    } 
    
    fit.AD <- optim(par = rep(0, ncol(Xcat)), logLik_AD, method = "BFGS") 
    se <- sqrt(diag(solve(hessian(func = function(x) logLik_AD(x)*sum(counts), x = fit.AD$par))))
    res <- list(est = fit.AD$par, se = se)
  }
  return(res)
}



