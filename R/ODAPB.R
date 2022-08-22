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


ODAPB.steps <- c('initialize','derive','estimate','synthesize')


# ODAPB code notes:
# 1.) Including options for Poisson, ZT-Poisson, quasi-Poisson, ZT-quasi-Poisson, Hurdle
#     - Lead site decides whether to use quasi-Poisson based on checking OD using their data. (can also return OD estimate in ODAPB.initialize if desired.)
#     - dist variable in each function specifies which Poisson dist to use.
# 2.) ODAPB methods include option for an offset variable. This is different from ODAL code, need to take this into account
#     (I think I did this correctly, but should be checked.)
# 3.) Will write separate functions for ODAH, since user can specify different sets of ipdata for
#     each component of hurdle model.


## write my own ztpoisson and hurdle to avoid import: countreg, as countreg is not on CRAN...
my.ztpoisson.loglik <- function(betas, X, Y, offset){
  design = as.matrix(X)
  betas <- as.matrix(betas)
  lp <- offset+c(design%*%betas)
  lambda <- exp(lp)
  sum((Y*lp - lambda - log(1-exp(-lambda)))*I(Y>0)) / length(Y[Y>0])   
}

my.ztpoisson.fit <- function(X, Y, offset){
  fn <- function(betas) - my.ztpoisson.loglik(betas, X, Y, offset)
  res <- optim(rep(0, ncol(X)), fn, hessian=T)
  return(list(b=res$par, b.var=diag(solve(res$hessian))/nrow(X)))
}



#' @useDynLib pda
#' @title ODAPB initialize
#' 
#' @usage ODAPB.initialize(ipdata,control,config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#' 
#' @references TBD
#' @return init
#' @keywords internal
ODAPB.initialize <- function(ipdata, control, config){
  # install.packages("countreg", repos="http://R-Forge.R-project.org")
  # dist <- control$dist
  family <- control$family
  # if(!("time" %in% colnames(ipdata))) {
  #   ipdata$time <- 1
  # }
  outcome <- ipdata$outcome
  offset <- ipdata$offset
  X <- as.matrix(ipdata[,-c(1,2)]) 
  L_gradient_meat = function(betas,X,Y,offset) {
    design <- as.matrix(X)
    betas = as.matrix(betas)
    Z = exp(design%*%betas)
    t(offset + c(Y-Z)*design) %*% (offset + c(Y-Z)*design)
  }
  #second-order gradient
  Lgradient2 <- function(betas, X, offset){
    design <- as.matrix(X)
    betas <- as.matrix(betas)
    t(- exp(offset + c(design%*%betas))*design)%*%design/nrow(X)
  }
  if (family == "poisson") {
    fit_i <- glm(outcome ~ 0+. -offset, data = ipdata, family = "poisson"(link = "log"), offset = offset)
    betas = fit_i$coef
    mat_L1 = L_gradient_meat(betas,X,outcome,offset)
    derivatives_i <- pdaGet(paste0(config$site_id,'_derive'),config)
    nlocal = derivatives_i$site_size
    
    mat_L2 = Lgradient2(betas,X,offset)*nlocal
    inv_L2 = solve.default(mat_L2)
    
    
    N <- 0
    for(site_i in control$sites){
      derivatives_i <- pdaGet(paste0(site_i,'_derive'),config)
      N <- N + derivatives_i$site_size
    }
    out = diag(inv_L2%*%mat_L1%*%inv_L2*(nlocal/N))
    init <- list(site = config$site_id,
                 site_size = nrow(ipdata),
                 bhat_i = fit_i$coef,
                 Vhat_i =  diag(vcov(fit_i)), 
                 phihat_i = 1)  
  }
  return(init)
}


#' @useDynLib pda
#' @title ODAPB derivatives
#' 
#' @usage ODAPB.derive(ipdata,control,config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#'
#' @return derivatives list(site = config$site_id, site_size = nrow(ipdata), logL_D1 = logL_D1, logL_D2 = logL_D2) 
#' @keywords internal
ODAPB.derive <- function(ipdata, control, config){
  family = control$family
  # dist <- control$dist
  # if(!("time" %in% colnames(ipdata))) {
  #   ipdata$time <- 1
  # }
  
  if (family == "poisson") {
    # data sanity check ...
    px <- ncol(ipdata) - 2  # number of covariates incl. intercept
    bhat <- rep(0, px)
    vbhat <- rep(0, px)     
    for(site_i in control$sites){
      init_i <- pdaGet(paste0(site_i,'_initialize'),config)
      bhat <- rbind(bhat, init_i$bhat_i)
      vbhat <- rbind(vbhat, init_i$Vhat_i)
    }
    bhat <- bhat[-1,] 
    vbhat <- vbhat[-1,]
    betameta = apply(bhat/vbhat,2,function(x){sum(x, na.rm = T)})/apply(1/vbhat,2,function(x){sum(x, na.rm = T)})
    vmeta = 1/apply(1/vbhat,2,function(x){sum(x, na.rm = T)})
    
    bbar <- betameta 
    
    # 1st and 2nd derivatives
    outcome <- ipdata$outcome
    offset <- ipdata$offset
    X <- as.matrix(ipdata[,-c(1,2)]) 
    # Getting rid of first two columns in ipdata (outcome and time) to define covariate matrix
    
    expit = function(x){1/(1 + exp(-x))}
    
    #first order gradient
    Lgradient <- function(betas, X, Y, offset){
      design <- as.matrix(X)
      betas <- as.matrix(betas)
      t(Y - exp(offset+design%*%betas)) %*% design / length(Y)
    }
    
    #second-order gradient
    Lgradient2 <- function(betas, X, offset){
      design <- as.matrix(X)
      betas <- as.matrix(betas)
      t(- exp(offset + c(design%*%betas))*design)%*%design/nrow(X)
    }
    
    logL_D1 <- Lgradient(bbar, X, outcome, offset)
    logL_D2 <- Lgradient2(bbar, X, offset)
    derivatives <- list(
      site = config$site_id, 
      site_size = nrow(ipdata),
      logL_D1 = logL_D1,
      logL_D2 = logL_D2)
  }
  
  return(derivatives)
}


#' @useDynLib pda
#' @title PDA surrogate estimation
#' 
#' @usage ODAPB.estimate(ipdata,control,config)
#' @param ipdata local data in data frame (generated in \code{pda})
#' @param control PDA control
#' @param config cloud configuration
#' 
#' @details  construct and solve surrogate logL at the master/lead site
#' @return  list(btilde = sol$par, Htilde = sol$hessian, site=control$mysite, site_size=nrow(ipdata))
#' @keywords internal
ODAPB.estimate <- function(ipdata,control,config) {
  family <- control$family
  # dist <- control$dist
  # if(!("time" %in% colnames(ipdata))) {
  #   ipdata$time <- 1
  # }
  
  if (family == "poisson") {
    # data sanity check ...
    Y <- ipdata$outcome
    outcome <- Y
    offset <- ipdata$offset
    X <- as.matrix(ipdata[,-c(1,2)]) 
    # Getting rid of first two columns in ipdata (outcome and time) to define covariate matrix
    px <- ncol(X)
    
    ######################################################
    #likelihood function for Poisson regression.
    Lik <- function(betas, X, Y, offset){
      design <- as.matrix(X)
      betas <- as.matrix(betas)
      lambda <- offset+c(design%*%betas)
      sum(Y*lambda - exp(lambda))/length(Y)
    }
    
    expit <- function(x){1/(1+exp(-x))}
    
    # download derivatives of other sites from the cloud
    # calculate 2nd order approx of the total logL
    logL_all_D1 <- rep(0, px)
    logL_all_D2 <- matrix(0, px, px)
    N <- 0
    for(site_i in control$sites){
      derivatives_i <- pdaGet(paste0(site_i,'_derive'),config)
      logL_all_D1 <- logL_all_D1 + derivatives_i$logL_D1*derivatives_i$site_size
      logL_all_D2 <- logL_all_D2 + derivatives_i$logL_D2*derivatives_i$site_size
      N <- N + derivatives_i$site_size
    }
    
    # initial beta
    bbar <- control$beta_init  # derivatives_i$b_meta
    
    #first order gradient
    Lgradient <- function(betas, X, Y, offset){
      design <- as.matrix(X)
      betas <- as.matrix(betas)
      t(Y - exp(offset+c(design%*%betas)))%*%design/length(Y)
    }
    
    #second-order gradient
    Lgradient2 <- function(betas, X, offset){
      design <- as.matrix(X)
      betas <- as.matrix(betas)
      t(-exp(offset+c(design%*%betas))*design)%*%design/nrow(X)
    }
    
    #second-order surogate likelihood, suppose the local data are stored in Xlocal, Ylocal
    # Y <- ipdata$outcome
    # t <- ipdata$time
    n1 <- length(Y)
    logL_tilde <- function(beta){
      - (Lik(beta,X, outcome, offset) + (logL_all_D1/N - Lgradient(bbar, X, outcome, offset))%*%beta+
           t(beta-bbar)%*%(logL_all_D2/N - Lgradient2(bbar, X, offset))%*%(beta-bbar) / 2)
    }
    
    # optimize the surrogate logL
    sol <- optim(par = bbar,
                 fn = logL_tilde,
                 # gr = logL_tilde_D1,
                 hessian = TRUE,
                 control = list(maxit=control$optim_maxit))
    
    if (family == "poisson") {
      surr <- list(btilde = sol$par, Htilde = sol$hessian,
                   site=config$site_id, site_size=nrow(ipdata))
    }
    ######################################################
  }
  
  return(surr)
}



ODAPB.synthesize <- function(ipdata, control, config) {
  family <- control$family
  
  px <- length(control$risk_factor)
  K <- length(control$sites)
  btilde_wt_sum <- rep(0, px)
  wt_sum <- rep(0, px)  
  
  for(site_i in control$sites){
    surr_i <- pdaGet(paste0(site_i,'_estimate'),config)
    btilde_wt_sum <- btilde_wt_sum + surr_i$Htilde %*% surr_i$btilde
    wt_sum <- wt_sum + surr_i$Htilde
  }
  
  # inv-Var weighted average est, and final Var = average Var-tilde
  btilde <- solve(wt_sum, btilde_wt_sum)
  Vtilde <- solve(wt_sum) * K
  
  message("all surrogate estimates synthesized, no need to broadcast! ")
  return(list(btilde = btilde,
              Vtilde = Vtilde))
  
}

