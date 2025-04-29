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

# Rcpp::sourceCpp('pda/src/rcpp_coxph.cpp')
# ODAC.steps<-c('initialize','deriveUWZ','derive','estimate','synthesize')
# ODAC.steps<-c('initialize','derive', 'estimate','synthesize')
# ODAC.family<-'cox'



#' @useDynLib pda
#' @title ODAC initialize
#' 
#' @usage ODAC.initialize(ipdata, control, config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#' 
#' @references  Rui Duan, et al. "Learning from local to global: An efficient distributed algorithm for modeling time-to-event data". 
#'               Journal of the American Medical Informatics Association, 2020, https://doi.org/10.1093/jamia/ocaa044
#'              Chongliang Luo, et al. ODACH: a one-shot distributed algorithm for Cox model with heterogeneous multi-center data. Sci Rep. 2022 Apr 22;12(1):6627. \doi{10.1038/s41598-022-09069-0}. 
#' @return  list(T_i = T_i, bhat_i = fit_i$coef, Vhat_i = summary(fit_i)$coef[,2]^2, site=control$mysite, site_size= nrow(ipdata))
#' @keywords internal
ODAC.initialize <- function(ipdata,control,config){
  if(control$heterogeneity == FALSE){
    T_i <- sort(unique(ipdata$time[ipdata$status==TRUE]))
  }else{
    T_i <- NA
  }
  
  # handle data degeneration (e.g. missing categories in some site). This could be in pda()?
  px = ncol(ipdata) - 2
  col_deg = apply(ipdata[,-c(1:2)],2,var)==0    # degenerated X columns...
  ipdata_i = ipdata[,-(which(col_deg)+2),with=F]
  
  # formula_i <- as.formula(paste("Surv(time, status) ~", paste(control$risk_factor[!col_deg], collapse = "+")) ) 
  fit_i <- tryCatch(survival::coxph(survival::Surv(time, status) ~ ., data=ipdata), error=function(e) NULL)
  
  if(!is.null(fit_i)){
    # for degenerated X, coef=0, var=Inf
    bhat_i = rep(0,px)
    Vhat_i = rep(Inf,px) 
    bhat_i[!col_deg] <- fit_i$coef
    Vhat_i[!col_deg] <- summary(fit_i)$coef[,"se(coef)"]^2 # not as glm, coxph summary can keep NA's! but vcov fills 0's!  
    
    init <- list(T_i = T_i,
                 bhat_i = bhat_i,
                 Vhat_i = Vhat_i,   #   
                 site = config$site_id,
                 site_size = nrow(ipdata))
    init$Vhat_i[init$Vhat_i==0] = NA
  } else{
    init <- list(T_i = T_i,
                 bhat_i = NA,
                 Vhat_i = NA,   
                 site = config$site_id,
                 site_size = nrow(ipdata))
  }
  
  return(init)
}



#' @useDynLib pda
#' @title Generate pda UWZ summary statistics before calculating derivatives
#' 
#' @usage ODAC.deriveUWZ(ipdata, control, config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#' @import Rcpp  
#' 
#' @return  list(T_all=T_all, b_meta=b_meta, site=control$mysite, site_size = nrow(ipdata), U=U, W=W, Z=Z, logL_D1=logL_D1, logL_D2=logL_D2)
#' @keywords internal
ODAC.deriveUWZ <- function(ipdata,control,config) {
  px <- ncol(ipdata) - 2
  # decide if doing ODAC derivatives 1st substep (calculate summary stats U, W, Z) 
  # or 2nd substep (calculate derivatives logL_D1, logL_D2)
  # collect event time pts and meta est from the cloud
  T_all <- c()
  # bhat_wt_sum <- rep(0, px)
  # wt_sum <- rep(0, px)     # cov matrix?
  for(site_i in control$sites){
    init_i <- pdaGet(paste0(site_i,'_initialize'),config)
    T_all <- c(T_all, as.numeric(init_i$T_i)) # in case no event T
    #   if(!any(is.na(init_i$bhat_i))){
    #     # NA may occur is some site has some degenerated X, meta will not use these sites but later surr lik can use them
    #     bhat_wt_sum <- bhat_wt_sum + init_i$bhat_i / init_i$Vhat_i
    #     wt_sum <- wt_sum + 1 / init_i$Vhat_i  # cov matrix?
    #   }
  }
  # 
  T_all <- sort(unique(T_all))
  nt <- length(T_all)
  # b_meta <- bhat_wt_sum / wt_sum
  # bbar <- b_meta
  ## (bug fix: ChongliangLuo, 20250105) 
  ## b_meta has already calculated from pdaSync lines 856-869 and stored in control
  bbar = control$beta_init
  
  # add fake data points to help calculate the summary stats in risk sets ar each time pts
  # t_max <- max(ipdata$time)+1
  t_max <- max(c(ipdata$time, T_all))+1 # 20250105
  # generate dataframe in format expected by ODAC
  pfdata <- cbind(T_all, 0, matrix(0, nt, px))
  pfdata <- rbind(ipdata, pfdata, use.names=FALSE)
  pfdata <- pfdata[, 'interval':=cut(pfdata$time, breaks = c(T_all, t_max), labels = 1:nt, right=FALSE)][order(pfdata$interval),]
  ## 20250105 bug fix: a surv T before T_all[1] should not be counted in any risk set!!!!
  # pfdata$interval[is.na(pfdata$interval)] <- nt
  pfdata = pfdata[!is.na(pfdata$interval),]
  # X <- as.matrix(pfdata[, control$variables, with=F])
  X <- as.matrix(pfdata[,-c(1,2)][,-'interval']) 
  
  # summary stats: U, W, Z
  eXb <- c(exp(X %*% bbar))
  X2 <- X[,1]*X
  for(ix in 2:ncol(X)) X2 <- cbind(X2, X[,ix]*X)
  UWZ <- eXb * cbind(1, X, X2)
  # rcpp_aggregate() is a function written in rcpp for calculating column-wise (reverse) cumsum
  # credit to Dr Wenjie Wang
  UWZ <- rcpp_aggregate(x = UWZ, indices = pfdata$interval, cumulative = T, reversely = T)
  
  # since fake X=0, cumulative W and Z will be the same, 
  # but exp(Xb)=1, so need to remove cumulated ones from each time pts
  # print(nrow(UWZ) )
  # print(nt)
  U <- UWZ[,1] - c(nt:1)
  W <- UWZ[,2:(px+1)]
  Z <- array(UWZ[,-c(1:(px+1))], c(nt,px,px))
  
  # summary_stat
  derivatives <- list(T_all=T_all, b_init=bbar, site=config$site_id, site_size=nrow(ipdata), U=U, W=W, Z=Z)
  
  # broadcast to the cloud?
  return(derivatives)
  
}


#' @useDynLib pda
#' @title Generate pda UWZ derivatives
#' 
#' @usage ODAC.derive(ipdata, control, config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#' 
#' @details Calculate and broadcast 1st and 2nd order derivative at initial bbar
#'        for ODAC, this requires 2 substeps: 1st calculate summary stats (U, W, Z), 
#'        2nd calculate derivatives (logL_D1, logL_D2)
#'
#' @import Rcpp  
#' @return  list(T_all=T_all, b_meta=b_meta, site=control$mysite, site_size = nrow(ipdata), U=U, W=W, Z=Z, logL_D1=logL_D1, logL_D2=logL_D2)
#' @keywords internal
ODAC.derive <- function(ipdata,control,config){
  px <- ncol(ipdata) - 2
  
  if (control$heterogeneity == F){
    # decide if doing ODAC derivatives 1st substep (calculate summary stats U, W, Z) 
    # or 2nd substep (calculate derivatives logL_D1, logL_D2)
    # read and add up (U W Z) from all sites from the cloud
    for(site_i in control$sites){
      sumstat_i <- pdaGet(paste0(site_i,'_deriveUWZ'),config)
      
      if(site_i == control$sites[1]){
        U <- sumstat_i$U
        W <- sumstat_i$W
        Z <- sumstat_i$Z
      }else{
        U <- U + sumstat_i$U
        W <- W + sumstat_i$W
        Z <- Z + sumstat_i$Z
      }
    }
    
    # number of events in ipdata at each event time pts in T_all
    T_all <- sumstat_i$T_all
    # bug: diff floating digits than T_all so need to round(time, 4) before running pda...
    d <- c(table(c(ipdata[ipdata$status==T,time], T_all)) - 1) 
    
    # 1st and 2nd derivatives
    # X <- as.matrix(ipdata[ipdata$status==TRUE, control$variables, with=F])
    X <- as.matrix(ipdata[ipdata$status==TRUE, -c(1,2)])
    
    logL_D1 <- apply(X, 2, sum) - apply(d * W / U, 2, sum, na.rm=T)
    W2 <- array(NA, c(dim(W), px))
    for(ii in 1:px) W2[,,ii] <- W[,ii] * W
    logL_D2 <- apply(d * (W2 - U*Z) / U^2, c(2, 3), sum, na.rm=T)  
    derivatives <- list(T_all=T_all, b_init=sumstat_i$b_init, U=U, W=W, Z=Z, 
                        site=config$site_id, site_size = nrow(ipdata),
                        logL_D1=logL_D1, logL_D2=logL_D2)
  } else { # ODACH
    time <- ipdata$time
    status <- ipdata$status
    X <- as.matrix(ipdata[,-c(1,2)])
    n <- length(time)
    # px <- ncol(X)
    
    ## get the initial values, beta_bar (i.e., bbar), broadcasted by sites
    # bhat_wt_sum <- rep(0, px)
    # wt_sum <- rep(0, px)     # cov matrix?
    # 
    # for(site_i in control$sites){
    #   init_i <- pdaGet(paste0(site_i,'_initialize'),config)
    #   bhat_wt_sum <- bhat_wt_sum + init_i$bhat_i / init_i$Vhat_i
    #   wt_sum <- wt_sum + 1 / init_i$Vhat_i  # cov matrix?
    # }
    # b_meta <- bhat_wt_sum / wt_sum
    # bbar <- b_meta
    ## (bug fix: AliFarnudi, 20240801) 
    ## b_meta has already calculated from pdaSync lines 804-835 and stored in control
    bbar = control$beta_init
     
    
    hasTies <- any(duplicated(ipdata$time)) 
    if(hasTies){
      # rcpp function is negative logL...
      logL_D1 <- -rcpp_coxph_logL_gradient_efron(bbar, time = time, event = status, z = X) # / n
      logL_D2 <- -matrix(rcpp_coxph_logL_hessian(bbar, time = time, event = status, z = X), px, px) # / n
    } else {
      logL_D1 <- -rcpp_coxph_logL_gradient(bbar, time = time, event = status, z = X) # / n
      logL_D2 <- -matrix(rcpp_coxph_logL_hessian(bbar, time = time, event = status, z = X), px, px) # / n
    }
    
    derivatives <- list(b_init=bbar, 
                        site=config$site_id, site_size = nrow(ipdata),
                        logL_D1=logL_D1, logL_D2=logL_D2)
  }
  
  return(derivatives)
}


#' @useDynLib pda
#' @title PDA surrogate estimation
#' 
#' @usage ODAC.estimate(ipdata, control, config)
#' @param ipdata local data in data frame
#' @param control pda control
#' @param config cloud config
#' @import data.table
#' 
#' @details step-4: construct and solve surrogate logL at the master/lead site
#' @import Rcpp  
#' @return  list(btilde = sol$par, Htilde = sol$hessian, site=control$mysite, site_size=nrow(ipdata))
#' @keywords internal
ODAC.estimate <- function(ipdata,control,config) {
  # data sanity check ...
  time <- ipdata$time
  status <- ipdata$status
  X <- as.matrix(ipdata[,-c(1,2)])
  n <- length(time)
  px <- ncol(X)
  hasTies <- any(duplicated(ipdata$time))
  
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
  # bbar <- derivatives_i$b_meta
  bbar <- control$beta_init
  
  # logL at local site
  if(hasTies){
    # rcpp function is negative logL...
    logL_local <- function(beta) -rcpp_coxph_logL_efron(beta, time = time, event = status, z = X) # / n
    logL_local_D1 <- function(beta) -rcpp_coxph_logL_gradient_efron(beta, time = time, event = status, z = X) # / n
    logL_local_D2 <- function(beta) -matrix(rcpp_coxph_logL_hessian(beta, time = time, event = status, z = X), px, px) # / n
  } else {
    logL_local <- function(beta) -rcpp_coxph_logL(beta, time = time, event = status, z = X)  # / n
    logL_local_D1 <- function(beta) -rcpp_coxph_logL_gradient(beta, time = time, event = status, z = X) # / n
    logL_local_D2 <- function(beta) -matrix(rcpp_coxph_logL_hessian(beta, time = time, event = status, z = X), px, px) # / n
  }
  
  # surrogate log-L and its gradient
  logL_diff_D1 <- logL_all_D1 / N - logL_local_D1(bbar) / n
  logL_diff_D2 <- logL_all_D2 / N - logL_local_D2(bbar) / n
  logL_tilde <- function(b) -(logL_local(b) / n + sum(b * logL_diff_D1) + 1/2 * t(b-bbar) %*% logL_diff_D2 %*% (b-bbar))
  # logL_tilde_D1 <- function(b) -(logL_local_D1(b) / n + logL_diff_D1 + logL_diff_D2 %*% (b-bbar))
  
  # optimize the surrogate logL 
  sol <- optim(par = bbar, 
               fn = logL_tilde,
               # gr = logL_tilde_D1,
               hessian = TRUE,
               method = control$optim_method,
               control = list(maxit=control$optim_maxit))
  
  # var estimate: by inv hessian 
  setilde = sqrt(diag(solve(sol$hessian))/N)
  
  surr <- list(btilde = sol$par, setilde=setilde, # Htilde = sol$hessian, 
               site=config$site_id, site_size=nrow(ipdata) )
  return(surr)
}



#' @useDynLib pda
#' @title PDA synthesize surrogate estimates from all sites, optional
#' 
#' @usage ODAC.synthesize(ipdata, control, config)
#' @param ipdata local data in data frame
#' @param control pda control
#' @param config cloud config
#' 
#' @details Optional step-4: synthesize all the surrogate est btilde_i from each site, if step-3 from all sites is broadcasted
#' @import Rcpp  
#' @return  list(btilde=btilde,  Vtilde=Vtilde)
#' @keywords internal
ODAC.synthesize <- function(ipdata,control,config) {
  px <- length(control$risk_factor)
  K <- length(control$sites)
  btilde = rep(0, px)
  setilde = rep(0, px) # cov matrix? 
  
  for(site_i in control$sites){
    surr_i <- pdaGet(paste0(site_i,'_estimate'),config)
    btilde = cbind(btilde, surr_i$btilde)
    setilde = cbind(setilde, surr_i$setilde) 
  }
  b_wt_sum <- rowSums(btilde / (setilde^2), na.rm=T)  
  wt_sum <- rowSums(1 / (setilde^2), na.rm=T) 
                         
  # inv-Var weighted average est, and final Var = average Var-tilde
  btilde <- b_wt_sum / wt_sum  
  setilde <- sqrt(K / wt_sum)  
  
  message("all surrogate estimates synthesized, no need to broadcast! ")
  return(list(btilde=btilde, setilde=setilde ))
}
