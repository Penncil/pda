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


ODAH.steps <- c('initialize','derive','estimate','synthesize')


# ODAP code notes:
# 1.) Including options for Poisson, ZT-Poisson, quasi-Poisson, ZT-quasi-Poisson, Hurdle
#     - Lead site decides whether to use quasi-Poisson based on checking OD using their data. (can also return OD estimate in ODAP.initialize if desired.)
#     - dist variable in each function specifies which Poisson dist to use.
# 2.) ODAP methods include option for an offset variable. This is different from ODAL code, need to take this into account
#     (I think I did this correctly, but should be checked.)
# 3.) Will write separate functions for ODAH, since user can specify different sets of ipdata for
#     each component of hurdle model.


## write my own hurdle to avoid import: countreg, as countreg is not on CRAN...
my.hurdle.fit <- function(X_count, X_zero, Y, offset){
  res.count <- my.ztpoisson.fit(X_count[Y>0,], Y[Y>0], offset[Y>0])
  res.zero <- glm(y~0+., data=cbind(y=ifelse(Y==0, 0, 1),X_zero), family='binomial')
  return(list(b.count=res.count$b,
              b.count.var=res.count$b.var,
              b.zero=res.zero$coef,
              b.zero.var=diag(vcov(res.zero))))
}




#' @useDynLib pda
#' @title ODAH initialize
#' 
#' @usage ODAH.initialize(ipdata,control,config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#' 
#' @references TBD
#' @return init
#' @keywords internal
ODAH.initialize <- function(ipdata, control, config){
  # install.packages("countreg", repos="http://R-Forge.R-project.org")
  # dist <- control$dist
  family <- control$family
  if(family!='hurdle'){
    warning("currently only support family='hurdle'" )
    family<-'hurdle'
  }
   
  # if (family == "hurdle") {
    # library(countreg)
    # outcome <- ipdata$outcome
    # X_count <- as.matrix(subset(ipdata, select = colnames(ipdata) %in% control_hurdle$variables_hurdle_count))
    # X_zero <- as.matrix(subset(ipdata, select = colnames(ipdata) %in% control_hurdle$variables_hurdle_zero))
    # fit_i <- hurdle(outcome ~ X_count + offset(log(ipdata$time)) | X_zero, dist = "poisson", zero.dist = "binomial")
    # fml <- as.formula(paste0(control$outcome, '~', paste0(control$variables_hurdle_count, collapse='+'), 
    #                          ifelse(is.character(control$offset), paste0('+offset(', control$offset, ')|'), '|'),
    #                          paste0(control$variables_hurdle_zero, collapse='+') ))
    # fit_i <-  countreg::hurdle(fml, data=ipdata$ipdata, dist = "poisson", zero.dist = "binomial")
    fit_i <- my.hurdle.fit(ipdata$X_count, ipdata$X_zero, ipdata$ipdata[,control$outcome], ipdata$offset)
    px_count <- length(control$variables_hurdle_count) + 1
    px_zero <- length(control$variables_hurdle_zero) + 1
    init <- list(site = config$site_id,
                 site_size = nrow(ipdata$ipdata),
                 bhat_zero_i = fit_i$b.zero,      # coefficients$zero,
                 bhat_count_i = fit_i$b.count,    # coefficients$count,
                 Vhat_zero_i = fit_i$b.zero.var,  # diag(vcov(fit_i))[-c(1:px_count)],
                 Vhat_count_i = fit_i$b.count.var # diag(vcov(fit_i))[1:px_count]
                ) 
  # }
  
  return(init)
}


#' @useDynLib pda
#' @title ODAH derivatives
#' 
#' @usage ODAH.derive(ipdata,control,config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#'
#' @return derivatives list(site = config$site_id, site_size = nrow(ipdata), 
#'                           logL_D1_zero = logL_D1_zero, logL_D1_count = logL_D1_count,
#'                           logL_D2_zero = logL_D2_zero, logL_D2_count = logL_D2_count) 
#' @keywords internal
ODAH.derive <- function(ipdata, control, config){
  family = control$family
  # dist <- control$dist
  # if(!("time" %in% colnames(ipdata))) {
  #   ipdata$time <- 1
  # }
  
  # if (family == "hurdle") {
    # X_count <- as.matrix(subset(ipdata, select = colnames(ipdata) %in% control$variables_hurdle_count))
    # X_zero <- as.matrix(subset(ipdata, select = colnames(ipdata) %in% control$variables_hurdle_zero))
    X_count <- ipdata$X_count
    X_zero <- ipdata$X_zero
    # data sanity check ...
    px_count <- ncol(X_count)  # number of covariates incl. intercept
    bhat_count <- rep(0, px_count)
    vbhat_count <- rep(0, px_count)
    px_zero <- ncol(X_zero)    # number of covariates incl. intercept
    bhat_zero <- rep(0, px_zero)
    vbhat_zero <- rep(0, px_zero)
    for(site_i in control$sites){
      init_i <- pdaGet(paste0(site_i,'_initialize'), config)
      bhat_zero <- rbind(bhat_zero, init_i$bhat_zero_i)
      vbhat_zero <- rbind(vbhat_zero, init_i$Vhat_zero_i)
      bhat_count <- rbind(bhat_count, init_i$bhat_count_i)
      vbhat_count <- rbind(vbhat_count, init_i$Vhat_count_i)
    }
    bhat_zero <- bhat_zero[-1,]
    vbhat_zero <- vbhat_zero[-1,]
    bhat_count <- bhat_count[-1,] 
    vbhat_count <- vbhat_count[-1,]
    
    betameta_zero <- apply(bhat_zero/vbhat_zero, 2, function(x){sum(x, na.rm = T)})/apply(1/vbhat_zero, 2, function(x){sum(x, na.rm = T)})
    betameta_count <- apply(bhat_count/vbhat_count, 2, function(x){sum(x, na.rm = T)})/apply(1/vbhat_count, 2, function(x){sum(x, na.rm = T)})
    vmeta_zero <- 1/apply(1/vbhat_zero, 2, function(x){sum(x, na.rm = T)})
    vmeta_count <- 1/apply(1/vbhat_count, 2, function(x){sum(x, na.rm = T)})
    
    bbar_zero <- betameta_zero
    bbar_count <- betameta_count
    
    # 1st and 2nd derivatives
    # outcome <- ipdata$ipdata$outcome
    outcome <- ipdata$ipdata[,control$outcome] 
    offset  <- ifelse(is.character(control$offset), ipdata$ipdata[,control$offset], 0)
    
    expit = function(x){1/(1 + exp(-x))}
    
    #first order gradients
    Lgradient_zero <- function(betas, X, Y){
      Y[Y>0] <- 1
      design <- as.matrix(X)
      betas <- as.matrix(betas)
      t(Y-expit(offset+c(design%*%betas)))%*%design/length(Y)
    }
    
    Lgradient_count <- function(betas, X, Y, offset){
      design <- as.matrix(X)
      betas <- as.matrix(betas)
      lambda <- exp(offset+c(design%*%betas))
      t(Y - lambda - lambda*exp(-lambda)/(1-exp(-lambda)))%*%(design*I(Y>0))/length(Y[Y>0])
    }
    
    #second-order gradients
    Lgradient2_zero <- function(betas, X){
      design <- as.matrix(X)
      betas <- as.matrix(betas)
      W <- expit(design%*%betas)
      t(c(-1*W*(1-W))*design)%*%design/nrow(X)
    }
    
    Lgradient2_count <- function(betas, X, Y, offset){
      design <- as.matrix(X)
      betas <- as.matrix(betas)
      lambda <- exp(offset+c(design%*%betas))
      t(c(-lambda - (exp(-lambda)*lambda*(1-exp(-lambda))*(1-lambda)+exp(-2*lambda)*exp(2*design%*%betas))/(1-exp(-lambda))^2)*design*I(Y>0))%*%(design*I(Y>0))/length(Y[Y>0])
    }
    
    logL_D1_zero <- Lgradient_zero(bbar_zero, X_zero, outcome)
    logL_D1_count <- Lgradient_count(bbar_count, X_count, outcome, offset)
    logL_D2_zero <- Lgradient2_zero(bbar_zero, X_zero)
    logL_D2_count <- Lgradient2_count(bbar_count, X_count, outcome, offset)
    
    derivatives <- list(
      site = config$site_id, 
      site_size = length(outcome), # nrow(ipdata),
      site_size_nonzero = sum(outcome>0),
      logL_D1_zero = logL_D1_zero,
      logL_D1_count = logL_D1_count,
      logL_D2_zero = logL_D2_zero,
      logL_D2_count = logL_D2_count) 
  # }
  return(derivatives)
}


#' @useDynLib pda
#' @title PDA surrogate estimation
#' 
#' @usage ODAH.estimate(ipdata,control,config)
#' @param ipdata local data in  a list(ipdata, X_count, X_zero) 
#' @param control PDA control
#' @param config cloud configuration
#' 
#' @details  construct and solve surrogate logL at the master/lead site
#' @return  list(btilde = sol$par, Htilde = sol$hessian, site=control$mysite, site_size=nrow(ipdata))
#' @keywords internal
ODAH.estimate <- function(ipdata,control,config) {
  family <- control$family
  # dist <- control$dist
  # if(!("time" %in% colnames(ipdata))) {
  #   ipdata$time <- 1
  # }
  
  # if (family == "hurdle") {
    # data sanity check ...
    outcome <- ipdata$ipdata[,control$outcome]
    offset  <- ifelse(is.character(control$offset), ipdata$ipdata[,control$offset], 0)
    X_count <- ipdata$X_count
    X_zero <- ipdata$X_zero
    # X_count <- as.matrix(subset(ipdata, select = colnames(ipdata) %in% control_hurdle$variables_hurdle_count))
    # X_zero <- as.matrix(subset(ipdata, select = colnames(ipdata) %in% control_hurdle$variables_hurdle_zero))
    
    px_count <- ncol(X_count)
    px_zero <- ncol(X_zero)
    
    ######################################################
    #likelihood function for ZT-Poisson regression.
    
    Lik_zero = function(betas, X, Y){
      Y[Y>0] <- 1
      design <- as.matrix(X)
      betas <- as.matrix(betas)
      lp <- offset + c(design%*%betas)
      sum(lp*Y - log(1+exp(lp)))/length(Y)
    }
    
    Lik_count = function(betas, X, Y, offset){
      # design = as.matrix(cbind(1,X))
      design = as.matrix(X)
      betas <- as.matrix(betas)
      lp <- offset+c(design%*%betas)
      lambda <- exp(lp)
      sum((Y*lp - lambda - log(1-exp(-lambda)))*I(Y>0))/length(Y[Y>0])  #  - lgamma(Y+1)
    }
    
    expit <- function(x){1/(1+exp(-x))}
    
    # download derivatives of other sites from the cloud
    # calculate 2nd order approx of the total logL
    logL_all_D1_zero <- rep(0, px_zero)
    logL_all_D1_count <- rep(0, px_count)
    logL_all_D2_zero <- matrix(0, px_zero, px_zero)
    logL_all_D2_count <- matrix(0, px_count, px_count)
    N <- 0
    N_nonzero <- 0 # need to get number of observations with counts > 0 at each site for ZT Poisson component of hurdle
    
    for(site_i in control$sites){
      derivatives_i <- pdaGet(paste0(site_i,'_derive'),config)
      logL_all_D1_zero <- logL_all_D1_zero + derivatives_i$logL_D1_zero*derivatives_i$site_size
      logL_all_D1_count <- logL_all_D1_count + derivatives_i$logL_D1_count*derivatives_i$site_size
      logL_all_D2_zero <- logL_all_D2_zero + derivatives_i$logL_D2_zero*derivatives_i$site_size
      logL_all_D2_count <- logL_all_D2_count + derivatives_i$logL_D2_count*derivatives_i$site_size
      N <- N + derivatives_i$site_size
      N_nonzero <- N_nonzero + derivatives_i$site_size_nonzero
    }
    
    # initial beta (needs to be coded into control_hurdle)
    bbar_zero <- control$beta_zero_init  
    bbar_count <- control$beta_count_init
    
    #first order gradients
    Lgradient_zero <- function(betas, X, Y){
      Y[Y>0] <- 1
      design <- as.matrix(X)
      betas <- as.matrix(betas)
      t(Y-expit(design%*%betas))%*%design/length(Y)
    }
    
    Lgradient_count <- function(betas, X, Y, offset){
      design <- as.matrix(X)
      betas <- as.matrix(betas)
      lambda <- exp(offset+c(design%*%betas))
      t(Y - lambda - lambda*exp(-lambda)/(1-exp(-lambda)))%*%(design*I(Y>0))/length(Y[Y>0])
    }
    
    #second-order gradients
    Lgradient2_zero <- function(betas, X){
      # design <- as.matrix(cbind(1,X))
      design <- as.matrix(X)
      betas <- as.matrix(betas)
      W <- expit(design%*%betas)
      t(c(-1*W*(1-W))*design)%*%design/nrow(X)
    }
    
    Lgradient2_count <- function(betas, X, Y, offset){
      # design <- as.matrix(cbind(1,X))
      design <- as.matrix(X)
      betas <- as.matrix(betas)
      lambda <- exp(offset+c(design%*%betas))
      t(c(-lambda - (exp(-lambda)*lambda*(1-exp(-lambda))*(1-lambda)+exp(-2*lambda)*lambda^2)/(1-exp(-lambda))^2)*design*I(Y>0))%*%(design*I(Y>0))/length(Y[Y>0])
    }
    
    #second-order surogate likelihood, logistic component: suppose the local data are stored in Xlocal, Ylocal
    # Y <- ipdata$outcome
    # Y[Y>0] <- 1
    # Y <- pmin(outcome, 1)
    n1 <- length(outcome)
    logL_tilde_zero <- function(beta){
      - (Lik_zero(beta, X_zero, outcome) + sum((logL_all_D1_zero/N - (Lgradient_zero(bbar_zero, X_zero, outcome)))*beta) +
           c(t(beta-bbar_zero)%*%(logL_all_D2_zero/N - Lgradient2_zero(bbar_zero, X_zero))%*%(beta-bbar_zero)) / 2)
    }
    
    #second-order surogate likelihood, ZT-Poisson component: suppose the local data are stored in Xlocal, Ylocal
    # Y <- ipdata$outcome
    # n1 <- length(Y[Y>0])
    logL_tilde_count <- function(beta){
      - (Lik_count(beta, X_count, outcome, offset) + sum((logL_all_D1_count/N_nonzero - (Lgradient_count(bbar_count, X_count, outcome, offset)))*beta)+
              c(t(beta-bbar_count)%*%(logL_all_D2_count/N_nonzero - Lgradient2_count(bbar_count, X_count, outcome, offset))%*%(beta-bbar_count)) / 2)
    }
    
    # optimize the surrogate logL
    sol_zero <- optim(par = bbar_zero,
                      fn = logL_tilde_zero,
                      # gr = logL_tilde_D1_zero,
                      hessian = TRUE,
                      control = list(maxit=control$optim_maxit))
    
    sol_count <- optim(par = bbar_count,
                       fn = logL_tilde_count,
                       # gr = logL_tilde_D1_count,
                       hessian = TRUE,
                       control = list(maxit=control$optim_maxit))
    
    
    surr <- list(btilde_zero = sol_zero$par, btilde_count = sol_count$par,
                 Htilde_zero = sol_zero$hessian, Htilde_count = sol_count$hessian,
                 site=config$site_id, site_size=n1 )  # nrow(ipdata)
    ######################################################
  # }
  return(surr)
}



ODAH.synthesize <- function(ipdata, control, config) {
  family <- control$family
  # dist <- control$dist
  # if(!("time" %in% colnames(ipdata))) {
  #   ipdata$time <- 1
  # }

  # if (control == control_hurdle) { # if using hurdle model
  # if(family=='hurdle'){
    px_zero <- length(control$variables_hurdle_zero) + 1
    px_count <- length(control$variables_hurdle_count) + 1
    K <- length(control$sites)
    btilde_zero_wt_sum <- rep(0, px_zero) 
    btilde_count_wt_sum <- rep(0, px_count) 
    wt_sum_zero <- rep(0, px_zero)
    wt_sum_count <- rep(0, px_count)
    
    for(site_i in control$sites){
      surr_i <- pdaGet(paste0(site_i,'_estimate'),config)
      btilde_zero_wt_sum <- btilde_zero_wt_sum + surr_i$Htilde_zero %*% surr_i$btilde_zero
      btilde_count_wt_sum <- btilde_count_wt_sum + surr_i$Htilde_count %*% surr_i$btilde_count
      wt_sum_zero <- wt_sum_zero + surr_i$Htilde_zero
      wt_sum_count <- wt_sum_count + surr_i$Htilde_count
    }
    
    # inv-Var weighted average est, and final Var = average Var-tilde
    btilde_zero <- solve(wt_sum_zero, btilde_zero_wt_sum)
    btilde_count <- solve(wt_sum_count, btilde_count_wt_sum)
    Vtilde_zero <- solve(wt_sum_zero) * K
    Vtilde_count <- solve(wt_sum_count) * K
    
    message("all surrogate estimates synthesized, no need to broadcast! ")
    return(list(btilde_zero = btilde_zero,
                btilde_count = btilde_count,
                Vtilde_zero = Vtilde_zero,
                Vtilde_count = Vtilde_count))
  
  # }
}

