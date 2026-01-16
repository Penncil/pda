# Copyright 2020 Penn Computing Inference Learning (PennCIL) lab
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

# ODACH_CC.steps<-c('initialize','derive', 'estimate','synthesize')
# ODACH_CC.family<-'cox'


#' @useDynLib pda
#' @title ODACH_CC initialize
#' 
#' @usage ODACH_CC.initialize(ipdata, control, config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#' 
#' @references  Chongliang Luo, et al. "ODACH: A One-shot Distributed Algorithm for Cox model with Heterogeneous Multi-center Data".
#'               medRxiv, 2021, https://doi.org/10.1101/2021.04.18.21255694
#' @return  list(bhat_i = fit_i$coef, Vhat_i = summary(fit_i)$coef[,2]^2, site=control$mysite, site_size= nrow(ipdata))
#' @keywords internal
ODACH_CC.initialize <- function(ipdata,control,config){
  # coxph with case-cohort design
  # ipdata columns: time_in, time_out, status, subcohort, strata_id, and covariates
  # full_cohort_size = control$full_cohort_size[control$sites==config$site_id]
  px = ncol(ipdata) - 5
  
  # handle data degeneration (e.g. missing categories in some site). This could be in pda()?
  col_deg = apply(ipdata[,-c(1:5)],2,var)==0    # degenerated X columns...
  ipdata_i = ipdata[,-(which(col_deg)+5),with=F]
  ipdata_i$ID = 1:nrow(ipdata_i) # for running coxph/cch... 
  
  # times <- sort(unique(c(ipdata$time_in, ipdata$time_out)))
  # precision <- min(diff(times)) / 2
  # ipdata_i[ipdata_i$subcohort == 0, "time_in"] <- ipdata_i[ipdata_i$subcohort == 0, "time_out"] - precision


  ## 3 ways to do local est: cch, coxph with a tweak of the formula, and cch_pooled
  # to avoid numerical error using cch() indicated by Ali, we use coxph with a tweak of the formula...
  # generally cch, coxph and cch_pooled will generate almost identical b and close var (for continuous X, coxph has smaller S.E. than the other two)
  # but coxph only works for Prentice wt, so will look into it later (and may revert to cch_pooled...)
  formula_i <- as.formula(paste("Surv(time_in, time_out, status) ~", 
                                paste(control$risk_factor[!col_deg], collapse = "+"), 
                                "+ strata(strata_id) + cluster(ID)"))
  fit_i <- tryCatch(survival::coxph(formula_i, data=ipdata_i, robust=T), error=function(e) NULL) 
  
  if(!is.null(fit_i)){
    # for degenerated X, coef=0, var=Inf
    bhat_i = rep(0,px)
    Vhat_i = rep(Inf,px) 
    bhat_i[!col_deg] <- fit_i$coef
    Vhat_i[!col_deg] <- summary(fit_i)$coef[,"se(coef)"]^2 # don't use robust var diag(fit_i$var)
    
    init <- list(bhat_i = bhat_i,
                 Vhat_i = Vhat_i,  
                 site = config$site_id,
                 site_size = nrow(ipdata),
                 # full_cohort_size = full_cohort_size, 
                 method = control$method)
    # init$Vhat_i[init$Vhat_i==0] = NA # 20250106
  } else{
    warning('survival::coxph() failed!!!')
    init <- list(bhat_i = rep(0,px),
                 Vhat_i = rep(Inf,px),  
                 S_i = NA,
                 site = config$site_id,
                 site_size = nrow(ipdata),
                 # full_cohort_size = full_cohort_size, 
                 method = control$method)
  }
  return(init)
}
 


get_logL_D1 <- function(coxph_fit){
  grad <- colSums(stats::residuals(coxph_fit, type = "score"))
  logL_D1 <- as.vector(grad)
  return(logL_D1)
}

get_logL_D2 <- function(coxph_fit){
  det <- survival::coxph.detail(coxph_fit, riskmat = FALSE)
  I_total <- apply(det$imat, c(1, 2), sum)
  logL_D2 <- -as.matrix(I_total)
  dimnames(logL_D2) <- NULL
  return(logL_D2)
}

#' @useDynLib pda
#' @title Generate pda derivatives
#' 
#' @usage ODACH_CC.derive(ipdata, control, config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#' 
#' @details Calculate and broadcast 1st and 2nd order derivative at initial bbar 
#'
#' @import Rcpp  
#' @return  list(bbar=bbar, site=control$mysite, site_size = nrow(ipdata), logL_D1=logL_D1, logL_D2=logL_D2)
#' @keywords internal
ODACH_CC.derive <- function(ipdata,control,config){
  px <- ncol(ipdata) - 5
  

  ## grad and hess
  bbar = control$beta_init
  
  # cc_prep = prepare_case_cohort(ipdata)
  # full_cohort_size = control$full_cohort_size[control$sites==config$site_id]


  # logL_D1 <- rep(0, px) # is it 0？
  # logL_D2 <- matrix(0, px, px)
  # logL_D1 <- rcpp_cc_grad_plk(beta = bbar,  
  #                             X = cc_prep$X,
  #                             failure_position = cc_prep$failure_position,
  #                             failure_num = cc_prep$failure_num,
  #                             risk_sets = cc_prep$risk_sets)
  # logL_D2 <- rcpp_cc_hess_plk(beta = bbar,  
  #                             X = cc_prep$X,
  #                             failure_num = cc_prep$failure_num,
  #                             risk_sets = cc_prep$risk_sets)
  
  # handle data degeneration (e.g. missing categories in some site). This could be in pda()?
  # col_deg = apply(ipdata[,-c(1:5)],2,var)==0    # degenerated X columns...
  # ipdata_i = ipdata[,-(which(col_deg)+5),with=F]
  # ipdata_i$ID = 1:nrow(ipdata_i) # for running coxph/cch...  
  
  # times <- sort(unique(c(ipdata$time_in, ipdata$time_out)))
  # precision <- min(diff(times)) / 2
  # ipdata_i[ipdata_i$subcohort == 0, "time_in"] <- ipdata_i[ipdata_i$subcohort == 0, "time_out"] - precision
  ## get intermediate (sandwich meat) for robust variance est of ODACH_CC  
  # formula_i <- as.formula(paste("Surv(time_in, time_out, status) ~", paste(control$risk_factor[!col_deg], collapse = "+"), '+ cluster(ID)')) 
  # fit_i <- tryCatch(survival::coxph(formula_i, data=ipdata_i, robust=T, init=bbar[!col_deg], iter=0), error=function(e) NULL) # 20250326: init/iter trick

  # times <- sort(unique(c(ipdata$time_in, ipdata$time_out)))
  # precision <- min(diff(times)) / 2
  # ipdata[ipdata$subcohort == 0, "time_in"] <- ipdata[ipdata$subcohort == 0, "time_out"] - precision
  ipdata$ID = 1:nrow(ipdata) # for running coxph/cch...  
  formula_i <- as.formula(
    paste("Surv(time_in, time_out, status) ~", 
    paste(control$risk_factor, collapse = "+"), '+ strata(strata_id) + cluster(ID)')) 
  fit_i <- survival::coxph(
    formula_i, 
    data=ipdata, 
    init=bbar, 
    method = "breslow",
    control = survival::coxph.control(iter.max = 0)
    )
  # y <- survival::Surv(ipdata$time_in, ipdata$time_out, ipdata$status)
  # X <- cc_prep$X
  # fit_i__ <- survival::coxph.fit(
  #   x = X,
  #   y = y, 
  #   strata=ipdata$strata,  
  #   init=bbar, 
  #   method = "breslow",
  #   rownames = as.character(seq_len(nrow(X))),
  #   control = survival::coxph.control(iter.max = 0),)
  # print(fit_i__$var)
  # print(logL_D1)

  # grad <- colSums(stats::residuals(fit_i, type = "score"))
  # logL_D1 <- as.vector(grad)
  logL_D1 <- get_logL_D1(fit_i)
  
  # print(as.vector(grad))
  # print(sum(abs(as.vector(grad)- logL_D1)))
  # print(as.vector(coef(fit_i)))
  # print(as.vector(coef(fit_i__)[!col_deg]))


  score_resid <- resid(fit_i, type = "score")  # n x p matrix  
  S_i <- matrix(0, px, px)   # this is the meat in sandwich var
  # S_i[!col_deg, !col_deg] <- crossprod(score_resid)
  S_i <- crossprod(score_resid)
  # print(S_i)
  # print(as.vector(crossprod(resid(fit_i, type = "score"))))
  # print(sum(abs(S_i- as.vector(crossprod(resid(fit_i__, type = "score"))))))
  # print(logL_D2)
  # print(-as.matrix(fit_i__$imat))
  # print(str(fit_i__))
  

  # det <- survival::coxph.detail(fit_i, riskmat = FALSE)
  # # print(dim(det$imat))
  # I_total <- apply(det$imat, c(1, 2), sum)
  # logL_D2 <- -as.matrix(I_total)
  # dimnames(logL_D2) <- NULL
  logL_D2 <- get_logL_D2(fit_i)
  
  # print(sum(abs(Hessian - logL_D2)))
  # I_surv <- Reduce(`+`, det$imat)
  # H_surv <- -I_surv
  # print(H_surv)

  # print(logL_D2[1])
  # print(as.matrix(solve(fit_i__$var))[1])
  # print(sum(abs(logL_D2[!col_deg, !col_deg]+as.matrix(solve(fit_i__$var[!col_deg, !col_deg])))))
  # stop("MES")
  # S_i[!col_deg, !col_deg] <- logL_D2[!col_deg, !col_deg] %*% fit_i$var %*% logL_D2[!col_deg, !col_deg] # Skhat in Yudong's note...
  
  derivatives <- list(bbar=bbar, 
                      site=config$site_id, 
                      site_size = nrow(ipdata), 
                      # full_cohort_size=full_cohort_size,
                      logL_D1=logL_D1, 
                      logL_D2=logL_D2,
                      S_i = S_i)
  
  return(derivatives)
}


#' @useDynLib pda
#' @title PDA surrogate estimation
#' 
#' @usage ODACH_CC.estimate(ipdata, control, config)
#' @param ipdata local data in data frame
#' @param control pda control
#' @param config cloud config
#' @import data.table
#' 
#' @details step-4: construct and solve surrogate logL at the master/lead site
#' @import Rcpp  
#' @return  list(btilde = sol$par, Htilde = sol$hessian, site=control$mysite, site_size=nrow(ipdata))
#' @keywords internal
ODACH_CC.estimate <- function(ipdata,control,config) {
  # data sanity check ... 
  px <- ncol(ipdata) - 5
  # hasTies <- any(duplicated(ipdata$time))
  
  # download derivatives of other sites from the cloud
  # calculate 2nd order approx of the total logL  
  logL_all_D1 <- rep(0, px)
  logL_all_D2 <- matrix(0, px, px)
  N <- 0
  for(site_i in control$sites){
    derivatives_i <- pdaGet(paste0(site_i,'_derive'),config)
    logL_all_D1 <- logL_all_D1 + derivatives_i$logL_D1
    logL_all_D2 <- logL_all_D2 + matrix(unlist(derivatives_i$logL_D2), px, px)
    N <- N + derivatives_i$site_size
  }
  
  # initial beta 
  bbar <- control$beta_init
  cc_prep = prepare_case_cohort(ipdata) 

  # logL at local site (mo negate or average)
  # logL_local <- function(beta) log_plk(beta, cc_prep)
  # logL_local_D1 <- function(beta) grad_plk(beta, cc_prep)
  # logL_local_D2 <- function(beta) hess_plk(beta, cc_prep)
  logL_local <- function(beta) rcpp_cc_log_plk(beta,  
                                               covariate = cc_prep$X,
                                               failure_position = cc_prep$failure_position,
                                               failure_num = cc_prep$failure_num,
                                               risk_sets = cc_prep$risk_sets)
  # logL_local_D1 <- function(beta) rcpp_cc_grad_plk(beta,  
  #                                                  X = cc_prep$X,
  #                                                  failure_position = cc_prep$failure_position,
  #                                                  failure_num = cc_prep$failure_num,
  #                                                  risk_sets = cc_prep$risk_sets)
  # logL_local_D2 <- function(beta) rcpp_cc_hess_plk(beta, 
  #                                                  X = cc_prep$X,
  #                                                  failure_num = cc_prep$failure_num,
  #                                                  risk_sets = cc_prep$risk_sets)
  
  # surrogate log-L and its gradient
  # logL_diff_D1 <- logL_all_D1 - logL_local_D1(bbar)  # / N / n (???? these loglik are not averaged)
  # logL_diff_D2 <- logL_all_D2 - logL_local_D2(bbar)  # / N / n
  logL_tilde <- function(b) -(logL_local(b) + sum(b * logL_diff_D1) + 1/2 * t(b-bbar) %*% logL_diff_D2 %*% (b-bbar)) #  / n
  # logL_tilde_D1 <- function(b) -(logL_local_D1(b) / n + logL_diff_D1 + logL_diff_D2 %*% (b-bbar)) 
  
  # times <- sort(unique(c(ipdata$time_in, ipdata$time_out)))
  # precision <- min(diff(times)) / 2
  # ipdata[ipdata$subcohort == 0, "time_in"] <- ipdata[ipdata$subcohort == 0, "time_out"] - precision
  ipdata$ID = 1:nrow(ipdata) # for running coxph/cch...  
  formula_i <- as.formula(
    paste("Surv(time_in, time_out, status) ~", 
    paste(control$risk_factor, collapse = "+"), '+ strata(strata_id) + cluster(ID)')) 
  fit_i <- survival::coxph(
    formula_i, 
    data=ipdata, 
    init=bbar, 
    method = "breslow",
    control = survival::coxph.control(iter.max = 0)
    )
  

  logL_local_D1_bbar <- get_logL_D1(fit_i)
  logL_diff_D1 <- logL_all_D1 - logL_local_D1_bbar

  logL_local_D2_bbar <- get_logL_D2(fit_i)
  logL_diff_D2 <- logL_all_D2 - logL_local_D2_bbar
  
  # col_deg = apply(ipdata[,-c(1:5)],2,var)==0    # degenerated X columns...
  # ipdata_i = ipdata[,-(which(col_deg)+5),with=F]
  # ipdata_i$ID = 1:nrow(ipdata_i) # for running coxph/cch... 
  # formula <- as.formula(paste("Surv(time_in, time_out, status) ~", 
  #                               paste(control$risk_factor[!col_deg], collapse = "+"), 
  #                               "+ strata(strata_id) + cluster(ID)"))
  # print(formula)
  
  # logL_tilde__ <- function(b){
  #   fit <- survival::coxph(
  #   formula_i, 
  #   data=ipdata, 
  #   init=b, 
  #   method = "breslow",
  #   control = survival::coxph.control(iter.max = 1)
  #   )
  #   res <-  -(fit$loglik[2] + sum(b * logL_diff_D1__) + 1/2 * t(b-bbar) %*% logL_diff_D2__ %*% (b-bbar))
  #   return(res)
  # }
#   logL_tilde__ <- function(b) {
#   fit <- try(
#     survival::coxph(
#       formula_i,
#       data = ipdata,
#       init = b,
#       method = "breslow",
#       control = survival::coxph.control(iter.max = 0)
#     ),
#     silent = TRUE
#   )
#   if (inherits(fit, "try-error")) {
#     return(1e100)
#   }
#   res <-  -(fit$loglik[2] + sum(b * logL_diff_D1__) + 1/2 * t(b-bbar) %*% logL_diff_D2__ %*% (b-bbar))
#   return(res)
# } 
  


  # optimize the surrogate logL 
  sol <- optim(par = bbar, 
               fn = logL_tilde,
               # gr = logL_tilde_D1,
               hessian = TRUE,
               method = control$optim_method,
               control = list(maxit=control$optim_maxit))
  
  # robust var estimate: see Yudong's note
  # setilde = sqrt(diag(solve(sol$hessian))/N)
  # hess of surrogate log-L at btilde, this is slightly diff than Yudong's, to avoid another iteration...
  # logL_tilde_D2 = logL_local_D2(bbar) + logL_diff_D2 
  logL_tilde_D2 = logL_local_D2_bbar + logL_diff_D2 
  # print(sum(abs(logL_tilde_D2-logL_tilde_D2__)))
  # put together
  # Stilde = solve(logL_tilde_D2) %*% control$S_i_sum %*% solve(logL_tilde_D2)
  solve_logL_tilde_D2 <- solve(logL_tilde_D2)
  Stilde = solve_logL_tilde_D2 %*% control$S_i_sum %*% solve_logL_tilde_D2
  setilde = sqrt(diag(Stilde))
  
  surr <- list(bbar=bbar, #full_cohort_size=full_cohort_size,  
               btilde = sol$par, setilde=setilde, Htilde = sol$hessian,
               site=config$site_id, site_size=nrow(ipdata))
  return(surr)
}



#' @useDynLib pda
#' @title PDA synthesize surrogate estimates from all sites, optional
#' 
#' @usage ODACH_CC.synthesize(ipdata, control, config)
#' @param ipdata local data in data frame
#' @param control pda control
#' @param config cloud config
#' 
#' @details Optional step-4: synthesize all the surrogate est btilde_i from each site, if step-3 from all sites is broadcasted
#' @import Rcpp  
#' @return  list(btilde=btilde,  Vtilde=Vtilde)
#' @keywords internal
ODACH_CC.synthesize <- function(ipdata,control,config) { 
  px <- length(control$risk_factor)
  K <- length(control$sites)
  btilde_wt_sum <- rep(0, px)
  wt_sum <- rep(0, px)     # cov matrix?
  
  for(site_i in control$sites){
    surr_i <- pdaGet(paste0(site_i,'_estimate'),config)
    btilde_wt_sum <- btilde_wt_sum + surr_i$Htilde %*% surr_i$btilde
    wt_sum <- wt_sum + surr_i$Htilde
  }
  
  # inv-Var weighted average est, and final Var = average Var-tilde
  btilde <- solve(wt_sum, btilde_wt_sum)
  Vtilde <- solve(wt_sum) * K
  
  message("all surrogate estimates synthesized, no need to broadcast! ")
  return(list(btilde=btilde, 
              Vtilde=Vtilde))
}
