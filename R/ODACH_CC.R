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
  ipdata$ID = 1:nrow(ipdata) # for running cch... 
  full_cohort_size = control$full_cohort_size[control$sites==config$site_id]
  formula <- as.formula(paste("Surv(time, status) ~", paste(control$variables, collapse = "+")))
  fit_i <- survival::cch(formula, data = ipdata, subcoh = ~subcohort, id = ~ID, 
                         cohort.size = full_cohort_size, method = control$method)
  # fit_i$var # summary(fit_i)$coef[,2]^2
  
  init <- list(bhat_i = fit_i$coef,
               Vhat_i = summary(fit_i)$coef[,"SE"]^2,   # not as glm, coxph summary can keep NA's! but vcov fills 0's!  
               site = config$site_id,
               site_size = nrow(ipdata),
               full_cohort_size = full_cohort_size, 
               method = control$method)
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
  # px <- ncol(ipdata) - 3
  
  bbar = control$beta_init
  full_cohort_size = control$full_cohort_size[control$sites==config$site_id]
  cc_prep = prepare_case_cohort(ipdata, control$method, full_cohort_size)
  # grad_plk() 
  logL_D1 <- grad_plk(bbar, cc_prep)
  # hess_plk()
  logL_D2 <- hess_plk(bbar, cc_prep)
  
  derivatives <- list(bbar=bbar, 
    site=config$site_id, site_size = nrow(ipdata), 
    full_cohort_size=full_cohort_size,
    logL_D1=logL_D1, logL_D2=logL_D2)
  
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
  cc_prep = prepare_case_cohort(ipdata, control$method, full_cohort_size) 
  
  # logL at local site
  logL_local <- function(beta) log_plk(beta, cc_prep)
  logL_local_D1 <- function(beta) grad_plk(beta, cc_prep)
  logL_local_D2 <- function(beta) hess_plk(beta, cc_prep)
  
  # surrogate log-L and its gradient
  logL_diff_D1 <- logL_all_D1 - logL_local_D1(bbar)  # / N / n
  logL_diff_D2 <- logL_all_D2 - logL_local_D2(bbar)  # / N / n
  logL_tilde <- function(b) -(logL_local(b) + sum(b * logL_diff_D1) + 1/2 * t(b-bbar) %*% logL_diff_D2 %*% (b-bbar)) #  / n
  # logL_tilde_D1 <- function(b) -(logL_local_D1(b) / n + logL_diff_D1 + logL_diff_D2 %*% (b-bbar))
  # cat(logL_diff_D1)
  # cat(logL_diff_D2)
  # cat(logL_tilde(bbar+0.2))
  
  # optimize the surrogate logL 
  sol <- optim(par = bbar, 
               fn = logL_tilde,
               # gr = logL_tilde_D1,
               hessian = TRUE,
               control = list(maxit=control$optim_maxit))
  
  surr <- list(bbar=bbar, full_cohort_size=full_cohort_size, 
               btilde = sol$par, Htilde = sol$hessian, site=config$site_id, site_size=nrow(ipdata))
  
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
