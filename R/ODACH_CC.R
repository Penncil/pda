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

# control <- list(project_name = 'IARC cancer study',
#                 step = 'initialize',
#                 sites = sites,
#                 heterogeneity = T,
#                 model = 'ODACH_CC', 
#                 method='Prentice', # 
#                 full_cohort_size=NA, # 
#                 family = 'cox',
#                 outcome = "status",
#                 variables = c('age', 'sex'),
#                 optim_maxit = 100,
#                 lead_site = 'site1',
#                 upload_date = as.character(Sys.time()) )

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
  full_cohort_size = control$full_cohort_size[control$sites==config$site_id]
  px = ncol(ipdata) - 3
  
  # handle data degeneration (e.g. missing categories in some site). This could be in pda()?
  col_deg = apply(ipdata[,-c(1:3)],2,var)==0    # degenerated X columns...
  ipdata_i = ipdata[,-(which(col_deg)+3),with=F]
  ipdata_i$ID = 1:nrow(ipdata_i) # for running coxph/cch... 
  
  # formula_i <- as.formula(paste("Surv(time, status) ~", paste(control$risk_factor[!col_deg], collapse = "+")))   
  # fit_i <- tryCatch(survival::cch(formula_i, data = ipdata_i, subcoh = ~subcohort, id = ~ID, 
  #                        cohort.size = full_cohort_size, method = control$method), error=function(e) NULL)  
  
  ## 3 ways to do local est: cch, coxph with a tweak of the formula, and cch_pooled
  # to avoid numerical error using cch() indicated by Ali, we use coxph with a tweak of the formula...
  # generally cch, coxph and cch_pooled will generate almost identical b and close var (for continuous X, coxph has smaller S.E. than the other two)
  # but coxph only works for Prentice wt, so will look into it later (and may revert to cch_pooled...)
  precision <- min(diff(sort(ipdata_i$time))) / 2 #  
  ipdata_i$time_in = 0
  ipdata_i[ipdata_i$subcohort == 0, "time_in"] <- ipdata_i[ipdata_i$subcohort == 0, "time"] - precision
  formula_i <- as.formula(paste("Surv(time_in, time, status) ~", paste(control$risk_factor[!col_deg], collapse = "+"), '+ cluster(ID)')) 
  fit_i <- tryCatch(coxph(formula_i, data=ipdata_i, robust=T), error=function(e) NULL) 
  
  if(!is.null(fit_i)){
    # for degenerated X, coef=0, var=Inf
    bhat_i = rep(0,px)
    Vhat_i = rep(Inf,px) 
    bhat_i[!col_deg] <- fit_i$coef
    Vhat_i[!col_deg] <- summary(fit_i)$coef[,"se(coef)"]^2 # dont's use robust var diag(fit_i$var)
    
    init <- list(bhat_i = bhat_i,
                 Vhat_i = Vhat_i,  
                 site = config$site_id,
                 site_size = nrow(ipdata),
                 full_cohort_size = full_cohort_size, 
                 method = control$method)
    # init$Vhat_i[init$Vhat_i==0] = NA # 20250106
  } else{
    warning('survival::coxph() failed!!!')
    init <- list(bhat_i = NA,
                 Vhat_i = NA,  
                 S_i = NA,
                 site = config$site_id,
                 site_size = nrow(ipdata),
                 full_cohort_size = full_cohort_size, 
                 method = control$method)
  }
  
  return(init)
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
  px <- ncol(ipdata) - 3
  
  # handle data degeneration (e.g. missing categories in some site). This could be in pda()?
  col_deg = apply(ipdata[,-c(1:3)],2,var)==0    # degenerated X columns...
  ipdata_i = ipdata[,-(which(col_deg)+3),with=F]
  ipdata_i$ID = 1:nrow(ipdata_i) # for running coxph/cch...  
  precision <- min(diff(sort(ipdata_i$time))) / 2 #  
  ipdata_i$time_in = 0
  ipdata_i[ipdata_i$subcohort == 0, "time_in"] <- ipdata_i[ipdata_i$subcohort == 0, "time"] - precision
  
  ## grad and hess
  bbar = control$beta_init
  full_cohort_size = control$full_cohort_size[control$sites==config$site_id]
  cc_prep = prepare_case_cohort(list(ipdata), control$method, full_cohort_size)
  # logL_D1 <- grad_plk(bbar, cc_prep)
  # logL_D2 <- hess_plk(bbar, cc_prep)
  logL_D1 <- rcpp_cc_grad_plk(beta = bbar, site_num = 1, 
               covariate_list = cc_prep$covariate_list,
               failure_position = cc_prep$failure_position,
               failure_num = cc_prep$failure_num,
               risk_sets = cc_prep$risk_sets,
               risk_set_weights = cc_prep$risk_set_weights)
  logL_D2 <- rcpp_cc_hess_plk(beta = bbar, site_num = 1,
                covariate_list = cc_prep$covariate_list,
                failure_position = cc_prep$failure_position,
                failure_num = cc_prep$failure_num,
                risk_sets = cc_prep$risk_sets,
                risk_set_weights = cc_prep$risk_set_weights)
  
  ## get intermediate (sandwich meat) for robust variance est of ODACH_CC  
  # fit_i <- tryCatch(coxph(formula_i, data=ipdata_i, robust=T), error=function(e) NULL) 
  formula_i <- as.formula(paste("Surv(time_in, time, status) ~", paste(control$risk_factor[!col_deg], collapse = "+"), '+ cluster(ID)')) 
  fit_i <- tryCatch(coxph(formula_i, data=ipdata_i, robust=T, init=bbar[!col_deg], iter=0), error=function(e) NULL) # 20250326: init/iter trick
  score_resid <- resid(fit_i, type = "score")  # n x p matrix  
  S_i = matrix(0, px, px)   # this is the meat in sandwich var
  S_i[!col_deg, !col_deg] <- crossprod(score_resid)
  # S_i[!col_deg, !col_deg] <- logL_D2[!col_deg, !col_deg] %*% fit_i$var %*% logL_D2[!col_deg, !col_deg] # Skhat in Yudong's note...
  
  derivatives <- list(bbar=bbar, 
    site=config$site_id, site_size = nrow(ipdata), 
    full_cohort_size=full_cohort_size,
    logL_D1=logL_D1, logL_D2=logL_D2,
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
  # time <- ipdata$time
  # status <- ipdata$status
  # X <- as.matrix(ipdata[,-c(1:3)])
  # n <- length(time)
  # px <- ncol(X)
  px <- ncol(ipdata) - 3
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
  # bbar <- derivatives_i$b_meta
  bbar <- control$beta_init
  full_cohort_size = control$full_cohort_size[control$sites==config$site_id]
  cc_prep = prepare_case_cohort(list(ipdata), control$method, full_cohort_size) 
  
  # logL at local site
  # logL_local <- function(beta) log_plk(beta, cc_prep)
  # logL_local_D1 <- function(beta) grad_plk(beta, cc_prep)
  # logL_local_D2 <- function(beta) hess_plk(beta, cc_prep)
  logL_local <- function(beta) rcpp_cc_log_plk(beta, site_num = 1, 
                             covariate_list = cc_prep$covariate_list,
                             failure_position = cc_prep$failure_position,
                             failure_num = cc_prep$failure_num,
                             risk_sets = cc_prep$risk_sets,
                             risk_set_weights = cc_prep$risk_set_weights)
  logL_local_D1 <- function(beta) rcpp_cc_grad_plk(beta, site_num = 1, 
                              covariate_list = cc_prep$covariate_list,
                              failure_position = cc_prep$failure_position,
                              failure_num = cc_prep$failure_num,
                              risk_sets = cc_prep$risk_sets,
                              risk_set_weights = cc_prep$risk_set_weights)
  logL_local_D2 <- function(beta) rcpp_cc_hess_plk(beta, site_num = 1,
                              covariate_list = cc_prep$covariate_list,
                              failure_position = cc_prep$failure_position,
                              failure_num = cc_prep$failure_num,
                              risk_sets = cc_prep$risk_sets,
                              risk_set_weights = cc_prep$risk_set_weights)
  
  # surrogate log-L and its gradient
  logL_diff_D1 <- logL_all_D1 - logL_local_D1(bbar)  # / N / n
  logL_diff_D2 <- logL_all_D2 - logL_local_D2(bbar)  # / N / n
  logL_tilde <- function(b) -(logL_local(b) + sum(b * logL_diff_D1) + 1/2 * t(b-bbar) %*% logL_diff_D2 %*% (b-bbar)) #  / n
  # logL_tilde_D1 <- function(b) -(logL_local_D1(b) / n + logL_diff_D1 + logL_diff_D2 %*% (b-bbar)) 
  
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
  logL_tilde_D2 = logL_local_D2(bbar) + logL_diff_D2 
  # put together
  Stilde = solve(logL_tilde_D2) %*% control$S_i_sum %*% solve(logL_tilde_D2)
  setilde = sqrt(diag(Stilde))
  
  surr <- list(bbar=bbar, full_cohort_size=full_cohort_size, 
               btilde = sol$par, setilde=setilde, Htilde = sol$hessian, site=config$site_id, site_size=nrow(ipdata))
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
