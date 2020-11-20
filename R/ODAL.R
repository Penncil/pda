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

# set in pda() ?
ODAL.steps <- c('initialize','derive','estimate','synthesize')
ODAL.family <- 'binomial'

#' @useDynLib pda
#' @title ODAL initialize
#' 
#' @usage ODAL.initialize(ipdata,control,config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#' 
#' @references Rui Duan, et al. "Learning from electronic health records across multiple sites: A communication-efficient and privacy-preserving distributed algorithm". 
#'                Journal of the American Medical Informatics Association, 2020, https://doi.org/10.1093/jamia/ocz199
#' @return init
#' @keywords internal
ODAL.initialize <- function(ipdata,control,config){
    fit_i <- glm(status ~ 0+., data=ipdata,family = "binomial"(link = "logit"))  
    init <- list(site = config$site_id,
                 site_size = nrow(ipdata),
                 bhat_i = fit_i$coef,
                 Vhat_i = diag(vcov(fit_i)))  # glm summary(fit_i)$coef[,2]^2 may omit NA's
  return(init)
}

#' @useDynLib pda
#' @title ODAL derivatives
#' 
#' @usage ODAL.derive(ipdata,control,config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#'
#' @return  list(site=config$site_id, site_size = nrow(ipdata), logL_D1=logL_D1, logL_D2=logL_D2)
#' @keywords internal
ODAL.derive <- function(ipdata,control,config){
  # data sanity check ...
    px <- ncol(ipdata) - 1  # X includes intercept
    # get b_meta as initial bbar
    bhat <- rep(0, px)
    vbhat <- rep(0, px)     # cov matrix?
    for(site_i in control$sites){
      init_i <- pdaGet(paste0(site_i,'_initialize'),config)
      bhat = rbind(bhat, init_i$bhat_i)
      vbhat = rbind(vbhat, init_i$Vhat_i)
    }
    bhat = bhat[-1,]
    vbhat = vbhat[-1,]
    
    #estimate from meta-analysis
    betameta = apply(bhat/vbhat,2,function(x){sum(x, na.rm = T)})/apply(1/vbhat,2,function(x){sum(x, na.rm = T)})
    vmeta = 1/apply(1/vbhat,2,function(x){sum(x, na.rm = T)})
    
    # b_meta <- betameta
    bbar <- betameta #b_meta
    
    # 1st and 2nd derivatives
    status <- ipdata$status
    X <- as.matrix(ipdata[,-1])

    expit = function(x){1/(1+exp(-x))}
    
    #first order gradient
    Lgradient = function(beta,X,Y){
      design = X  # cbind(1,X)
      t(Y-expit(design%*%t(t(beta))))%*%design/length(Y)
    }
    
    #second-order gradient
    Lgradient2 = function(beta,X){
      design = X # cbind(1,X)
      Z=expit(design%*%beta)
      t(c(-Z*(1-Z))*design)%*%design/nrow(X)
    }
    
    logL_D1 <- Lgradient(bbar,X,ipdata$status)
    logL_D2 <- Lgradient2(bbar,X)
    
    derivatives <- list(
      site=config$site_id, 
      site_size = nrow(ipdata),
      logL_D1=logL_D1,
      logL_D2=logL_D2)
  
  return(derivatives)
}


#' @useDynLib pda
#' @title PDA surrogate estimation
#' 
#' @usage ODAL.estimate(ipdata,control,config)
#' @param ipdata local data in data frame
#' @param control PDA control
#' @param config cloud configuration
#' 
#' @details step-3: construct and solve surrogate logL at the master/lead site
#' @return  list(btilde = sol$par, Htilde = sol$hessian, site=control$mysite, site_size=nrow(ipdata))
#' @keywords internal
ODAL.estimate <- function(ipdata,control,config) {
  # data sanity check ...
    status <- ipdata$status
    X <- as.matrix(ipdata[,-1])
    px <- ncol(X)
    
    ######################################################
    #likelihood function for logistic regression, the input X is a n*d matrix where
    #each patient has d covariates stored in each row.
    Lik = function(beta,X,Y){
      design = X # cbind(1,X)
      sum(Y*(design%*%t(t(beta)))-log(1+exp(design%*%t(t(beta)))))/length(Y)
    }
    
    expit = function(x){1/(1+exp(-x))}
    
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
      Lgradient = function(beta,X,Y){
        design = X  # cbind(1,X)
        t(Y-expit(design%*%t(t(beta))))%*%design/length(Y)
      }
      #second-order gradient
      Lgradient2 = function(beta,X){
        design = X # cbind(1,X)
        Z=expit(design%*%beta)
        t(c(-Z*(1-Z))*design)%*%design/nrow(X)
      }
      
      #first-order surogate likelihood, suppose the local data are stored in Xlocal, Ylocal
      Y = ipdata$status
      n1 = length(Y)
      logL_tilde = function(beta){
        - (Lik(beta,X, Y) + (logL_all_D1/N - Lgradient(bbar, X, Y))%*%beta+
             t(beta-bbar)%*%(logL_all_D2/N - Lgradient2(bbar, X))%*%(beta-bbar) / 2)
      }
      
      # optimize the surrogate logL
      sol <- optim(par = bbar,
                   fn = logL_tilde,
                   # gr = logL_tilde_D1,
                   hessian = TRUE,
                   control = list(maxit=control$optim_maxit))
      
    
    surr <- list(btilde = sol$par, Htilde = sol$hessian, site=config$site_id, site_size=nrow(ipdata))
    ######################################################
    
  return(surr)
}



#' @useDynLib pda
#' @title PDA synthesize surrogate estimates from all sites, optional
#' 
#' @usage ODAL.synthesize(ipdata,control,config)
#' @param ipdata local data in data frame
#' @param control pda control
#' @param config pda cloud configuration
#' 
#' @details Optional step-4: synthesize all the surrogate est btilde_i from each site, if step-3 from all sites is broadcasted
#'
#' @return  list(btilde=btilde,  Vtilde=Vtilde)
#' @keywords internal
ODAL.synthesize <- function(ipdata,control,config) {
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
