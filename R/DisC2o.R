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
DisC2o.steps <- c('PSinitialize','PSderive','PSestimate',
                  'OMinitialize','OMderive','OMestimate',
                  'AIPWestimate','synthesize')
# Rcpp::sourceCpp("pda/src/DisC2o.cpp")


#' @useDynLib pda
#' @title DisC2o PS initialize
#' 
#' @usage DisC2o.PSinitialize(ipdata,control,config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#' 
#' @references Tong J, et al. 2025. DisC2o-HD: Distributed causal inference with covariates shift for analyzing real-world high-dimensional data. Journal of Machine Learning Research. 2025;26(3):1-50.
#' @return init
#' @keywords internal
DisC2o.PSinitialize <- function(ipdata,control,config){
  # handle data degeneration (e.g. missing categories in some site). This could be in pda()?
  px = ncol(ipdata) - 2
  col_deg = apply(ipdata[, -c(1:3), with = FALSE], 2, var) == 0    # degenerated X columns...
  ipdata_i = ipdata[,-(which(col_deg)+3),with=F]
  
  # fit_i <- tryCatch(glm(status ~ 0+., data=ipdata_i, family = "binomial"(link = "logit")), error=function(e) NULL)
  # fit_i <- glm(status ~ 0+., data=ipdata,family = "binomial"(link = "logit"))  
  x_mat <- as.matrix(ipdata_i[, -c(1:3), with = FALSE])  # predictors
  y_vec <- ipdata_i[["treatment"]]                     # response as a vector
  
  fit0 <- cv.glmnet(x_mat, y_vec,  family = "binomial")
  lambda<-fit0$lambda.min
  fit_i <- glmnet(x_mat, y_vec,  family = "binomial",lambda=lambda)
  
  if(!is.null(fit_i)){
    # for degenerated X, coef=0, var=Inf
    bhat_i = rep(0,px)
    # Vhat_i = rep(Inf,px) 
    bhat_i[c(1,which(!col_deg)+1)] <- as.vector(coef(fit_i))
    #Vhat_i[!col_deg] <- diag(vcov(fit_i)) # summary(fit_i)$coef[,2]^2 may omit NA's
    init <- list(thetahat_i = bhat_i,
                 # Vhat_i = Vhat_i,
                 site = config$site_id,
                 site_size = nrow(ipdata))   
  } else{
    init <- list(thetahat_i = NA,
                 #Vhat_i = NA,   
                 site = config$site_id,
                 site_size = nrow(ipdata))
  }
  return(init)
}

#' @useDynLib pda
#' @title DisC2o_PS derivatives
#' 
#' @usage DisC2o.PSderive(ipdata,control,config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#'
#' @return  list(site=config$site_id, site_size = nrow(ipdata), logL_D1=logL_D1, logL_D2=logL_D2)
#' @keywords internal
DisC2o.PSderive <- function(ipdata,control,config){
  # data sanity check ...
    px <- ncol(ipdata) - 2  # X includes intercept 
    
    # get b_meta as initial bbar
    init_i <- pdaGet(paste0(control$lead_site,'_PSinitialize'),config)
    bbar <- init_i$thetahat_i 
    
    # 1st and 2nd derivatives
    treat <- ipdata$treatment
    X <- as.matrix(ipdata[,-c(1,2)])
 
    expit = function(x){1/(1+exp(-x))}
    
    #first order gradient
    GH_bio_ho <- function(beta, X,  Treat){
      n<-length(Treat)
      mean0 <-X%*%beta
      p <- expit(mean0)
      weight <- diag(c(p*(1-p)))
      I <- (t(X)%*%weight%*%X)/n
      S <- t(Treat-p)%*%X/n
      return(list(gradient=t(S), hessian=I))
    }
    
    fit<-GH_bio_ho(bbar, X, treat)
    
    logL_D1 <- fit$gradient
    logL_D2 <- fit$hessian
    
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
#' @usage DisC2o.PSestimate(ipdata,control,config)
#' @param ipdata local data in data frame
#' @param control PDA control
#' @param config cloud configuration
#' 
#' @return  list(btilde = sol$par, Htilde = sol$hessian, site=control$mysite, site_size=nrow(ipdata))
#' @keywords internal
DisC2o.PSestimate <- function(ipdata,control,config) {
  # data sanity check ...
    treat<- ipdata$treatment
    X <- as.matrix(ipdata[,-c(1,2)])
    px <- ncol(X)
    
    ######################################################
    #likelihood function for logistic regression, the input X is a n*d matrix where
    #each patient has d covariates stored in each row.
    Lik_bio <- function(beta, X, Treat){
      -mean(Treat * (X%*%beta)- log(1 + exp(X%*%beta)))
    }
    
    expit = function(x){1/(1+exp(-x))}
    
    #first order gradient
    GH_bio_ho <- function(beta, X,  Treat){
      n<-length(Treat)
      mean0 <-X%*%beta
      p <- expit(mean0)
      weight <- diag(c(p*(1-p)))
      I <- (t(X)%*%weight%*%X)/n
      S <- t(Treat-p)%*%X/n
      return(list(gradient=t(S), hessian=I))
    } 
    
    cv.adap_bio_ho<-function(betainit, 
                             betabar,  
                             X1,  
                             Treat, 
                             L1all, 
                             L2all=NULL,
                             nlambda=20,
                             nfolds=5,
                             method=c("CV","BIC")){
      n<-length(Treat)
      fit1<-GH_bio_ho(betabar, X1, Treat)
      fit<-GH_bio_ho(betainit, X1, Treat)
      if(is.null(L2all)){
        B <- fit$hessian
        atilde <--(fit$gradient+ L1all-fit1$gradient)-B%*%betainit
      }else{
        B <- fit$hessian+L2all-fit1$hessian
        atilde <- -(fit$gradient + L1all-fit1$gradient
                    -(L2all-fit1$hessian)%*%(betainit-betabar))-B%*%betainit
      }  
      
      lam.max <- max(abs((atilde + (B - diag(diag(B)))%*%betainit)/diag(B))[-1])
      lam.min <- 0.02*lam.max
      lam.seq <- exp(seq(log(lam.min), log(lam.max), length =nlambda ))
      
      if(method=="CV"){
        folds=cv.folds(n,nfolds)
        re<-NULL
        for (lambda in lam.seq) {
          se<-0
          for (k in 1:nfolds) {
            Treat.train=Treat[as.vector(unlist(folds[-k]))]
            X1.train=X1[as.vector(unlist(folds[-k])),]
            
            Treat.test=Treat[as.vector(unlist(folds[k]))]
            X1.test=X1[as.vector(unlist(folds[k])),]
            
            fit<-adap_bio_ho(betainit, betabar,  X1.train,  Treat.train,
                             L1all, L2all,lambda)
            beta<-fit$Beta_est_adap
            se=se+Lik_bio(beta, X1.test, Treat.test)
          }
          re<-c(re,se)
        }
        lambda<-lam.seq[which.min(re)]
        fit<-adap_bio_ho(betainit, betabar,  X1,  Treat, L1all, L2all,lambda)
        beta<-fit$Beta_est_adap
        msg<-fit$msg
      }
      
      if(method=="BIC"){
        re<-NULL
        for(lambda in lam.seq){
          fit <- adap_bio_ho(betainit, betabar,  X1,  Treat, L1all, L2all,lambda)
          beta <- fit$Beta_est_adap
          re <- c(re, 2*Lik_bio(beta, X1, Treat)+log(n)*sum(I(beta!=0)/n))
        }
        lambda <- lam.seq[which.min(re)]
        fit <- adap_bio_ho(betainit, betabar,  X1,  Treat,
                         L1all, L2all,lambda)
        beta<-fit$Beta_est_adap
        msg<-fit$msg
      }
      if(msg==0){
        cat("Successful converge!", "\n")
      }else{
        stop("Lambda is too small!", "\n")
      } 
      
      return(list(coeff=beta, msg=msg))  
    }
    
    
    #coordinate descent algorithm for one-shot
    adap_bio_ho<-function( betainit, 
                           betabar,  
                           X1,  
                           Treat, 
                           L1all, 
                           L2all=NULL,
                           lambda){
      fit1<-GH_bio_ho(betabar, X1, Treat)
      tol <- 1e-4
      out.loop <- 1
      out.Loop <- 100
      msg <- 0
      while (out.loop <= out.Loop) {
        fit<-GH_bio_ho(betainit, X1, Treat)
        if(is.null(L2all)){
          B <- fit$hessian
          atilde <--(fit$gradient+ L1all-fit1$gradient)-B%*%betainit
        }else{
          B <- fit$hessian+L2all-fit1$hessian
          atilde <- -(fit$gradient + L1all-fit1$gradient
                      -(L2all-fit1$hessian)%*%(betainit-betabar))-B%*%betainit
        }
        fitcd <- coordi.c(atilde, B, betainit, lambda, iter_max=100)
        beta.new <- fitcd$betainit
        if (fitcd$message==1){
          msg<-1
          break
        }
        dif <- sqrt(t(beta.new-betainit)%*%(beta.new-betainit))
        betainit <- beta.new
        if (abs(dif) < tol) {
          #cat("At OUTER loop", out.loop, "Successful converge", "\n")
          msg<-0
          break
        }
        if (out.loop == out.Loop) {
          #cat("OUTER loop maximum iteration reached", "\n")
          msg<-1
        }
        out.loop <- out.loop + 1
      }
      return(list(Beta_est_adap=beta.new, msg=msg))
    }
    
    
    #coordinate descent algorithm
    coordi.c<-function(atilde, B, betainit, lambda, iter_max=100){
      msg<-1
      fit<-coordi_cho(atilde, B, betainit, lambda,iter_max)
      if (fit$iter<= iter_max) {
        #cat("At INNER loop", loop, "Successful converge", "\n")
        msg <- 0
      } 
      return(list(betainit=fit$betainit, message=msg))
    } 
    
    # download derivatives of other sites from the cloud
    # calculate 2nd order approx of the total logL
      logL_all_D1 <- rep(0, px)
      logL_all_D2 <- matrix(0, px, px)
      N <- 0
      for(site_i in control$sites){
        derivatives_i <- pdaGet(paste0(site_i,'_PSderive'),config)
        logL_all_D1 <- logL_all_D1 + derivatives_i$logL_D1*derivatives_i$site_size
        logL_all_D2 <- logL_all_D2 + derivatives_i$logL_D2*derivatives_i$site_size
        N <- N + derivatives_i$site_size
      }
      
      logL_all_D1<-logL_all_D1/N
      logL_all_D2<-logL_all_D2/N
      
      # initial beta
      init_i <- pdaGet(paste0(control$lead_site,'_PSinitialize'),config)
      bbar <- init_i$thetahat_i 
      #bbar <- control$beta_init  # derivatives_i$b_meta
      
      fit<-cv.adap_bio_ho(bbar, 
                          bbar,  
                          X1=X,  
                          Treat=treat, 
                          logL_all_D1, 
                          logL_all_D2,
                          nlambda=20,
                          nfolds=5,
                          method="BIC")
      
    # Htilde = sol$hessian, 
    surr <- list(thetatilde = fit$coeff, 
                 #setilde=sqrt(diag(solve(sol$hessian))/N), 
                 site=config$site_id, 
                 site_size=nrow(ipdata))
    
  return(surr)
}


#' @useDynLib pda
#' @title DisC2o_OM initialize
#' 
#' @usage DisC2o.OMinitialize(ipdata,control,config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#' 
#' @return init
#' @keywords internal
DisC2o.OMinitialize <- function(ipdata,control,config){
  # handle data degeneration (e.g. missing categories in some site). This could be in pda()?
  px = ncol(ipdata) - 2
  col_deg = apply(ipdata[, -c(1:3), with = FALSE], 2, var) == 0    # degenerated X columns...
  ipdata_i = ipdata[,-(which(col_deg)+3),with=F]
  
  # fit_i <- tryCatch(glm(status ~ 0+., data=ipdata_i, family = "binomial"(link = "logit")), error=function(e) NULL)
  # fit_i <- glm(status ~ 0+., data=ipdata,family = "binomial"(link = "logit"))  
  x_mat <- as.matrix(ipdata_i[, -c(1:3), with = FALSE])  # predictors
  y_vec <- ipdata_i[["status"]]                     # response as a vector
  treat<-ipdata_i[["treatment"]]
  
  # function for computing IPW weights
  weight_treat_ho<-function(X, 
                            thetahat, 
                            Treat){
    mean1<-X%*%thetahat
    weight<-sqrt(Treat/exp(mean1))
    return(as.numeric(weight))
  }
  
  init_i <- pdaGet(paste0(control$lead_site,'_PSestimate'),config)
  thetahat<-init_i$thetatilde
  ipw_weight<-weight_treat_ho(cbind(1,x_mat), thetahat, treat)
  
  
  fit0 <- cv.glmnet(x_mat, y_vec, weights = ipw_weight, family = "gaussian")
  lambda<-fit0$lambda.min
  fit_i <- glmnet(x_mat, y_vec, lambda=lambda, weights = ipw_weight, family = "gaussian")
  
  
  if(!is.null(fit_i)){
    # for degenerated X, coef=0, var=Inf
    bhat_i = rep(0,px)
    # Vhat_i = rep(Inf,px) 
    bhat_i[c(1,which(!col_deg)+1)] <- as.vector(coef(fit_i))
    #Vhat_i[!col_deg] <- diag(vcov(fit_i)) # summary(fit_i)$coef[,2]^2 may omit NA's
    init <- list(bhat_i = bhat_i,
                 # Vhat_i = Vhat_i,
                 site = config$site_id,
                 site_size = nrow(ipdata))   
  } else{
    init <- list(bhat_i = NA,
                 #Vhat_i = NA,   
                 site = config$site_id,
                 site_size = nrow(ipdata))
  }
  return(init)
}

#' @useDynLib pda
#' @title DisC2o_OM derivatives
#' 
#' @usage DisC2o.OMderive(ipdata,control,config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#'
#' @return  list(site=config$site_id, site_size = nrow(ipdata), logL_D1=logL_D1, logL_D2=logL_D2)
#' @keywords internal
DisC2o.OMderive <- function(ipdata,control,config){
  # data sanity check ...
  px <- ncol(ipdata) - 1  # X includes intercept 
  
  # get b_meta as initial bbar
  init_i <- pdaGet(paste0(control$lead_site,'_OMinitialize'),config)
  bbar <- init_i$bhat_i 
  #bbar <- control$beta_init 
  
  # 1st and 2nd derivatives
  status <- ipdata$status
  X <- as.matrix(ipdata[,-c(1,2)])
  treat<-ipdata$treatment
  
  # function for computing IPW weights
  weight_treat_ho<-function(X, 
                            thetahat, 
                            Treat){
    mean1<-X%*%thetahat
    weight<-sqrt(Treat/exp(mean1))
    return(as.numeric(weight))
  }
  
  init_i <- pdaGet(paste0(control$lead_site,'_PSestimate'),config)
  thetahat<-init_i$thetatilde
  ipw_weight<-weight_treat_ho(X, thetahat, treat)
  
  
  #first and second order gradient
  GH_ls_ho <- function(beta,
                       X, 
                       weight, 
                       Y){
    Y<-weight*Y
    n<-length(Y)
    
    I <- (t(X)%*%X)/n
    S <- t(X)%*%(Y-X%*%beta)/n
    
    return(list(gradient=S, hessian=I))
  }
  
  
  fit<-GH_ls_ho(beta=bbar,
                X=X,
                weight=ipw_weight,
                Y=status)
  
  logL_D1 <- fit$gradient
  logL_D2 <- fit$hessian
  
  derivatives <- list(
    site=config$site_id, 
    site_size = nrow(ipdata),
    logL_D1=logL_D1,
    logL_D2=logL_D2)
  
  return(derivatives)
}


#' @useDynLib pda
#' @title DisC2o outcome model surrogate estimation
#' 
#' @usage DisC2o.OMestimate(ipdata,control,config)
#' @param ipdata local data in data frame
#' @param control PDA control
#' @param config cloud configuration
#' 
#' @return  list(btilde = sol$par, Htilde = sol$hessian, site=control$mysite, site_size=nrow(ipdata))
#' @keywords internal
DisC2o.OMestimate <- function(ipdata,control,config) {
  # data sanity check ...
  status <- ipdata$status
  X <- as.matrix(ipdata[,-c(1,2)])
  treat <- ipdata$treatment
  px <- ncol(X)
  
  ######################################################
  #likelihood function for logistic regression, the input X is a n*d matrix where
  #each patient has d covariates stored in each row.
  
  # function for computing IPW weights
  weight_treat_ho<-function(X, 
                            thetahat, 
                            Treat){
    mean1<-X%*%thetahat
    weight<-sqrt(Treat/exp(mean1))
    return(as.numeric(weight))
  }
  
  init_i <- pdaGet(paste0(control$lead_site,'_PSestimate'),config)
  thetahat <- init_i$thetatilde
  #thetahat<-control$thetahat
  ipw_weight<-weight_treat_ho(X, thetahat, treat)
  
  #first order gradient
  GH_ls_ho <- function(beta,
                       X, 
                       weight, 
                       Y){
    Y<-weight*Y
    n<-length(Y)
    
    I <- (t(X)%*%X)/n
    S <- t(X)%*%(Y-X%*%beta)/n
    
    return(list(gradient=S, hessian=I))
  }
  
  cv.adap_ls_ho<-function(betainit, 
                          betabar,  
                          X1,  
                          Y, 
                          L1all, 
                          L2all=NULL,
                          nlambda=20,
                          weight,
                          nfolds=5,
                          method=c("CV","BIC")){
    n<-length(Y)
    fit1<-GH_ls_ho(betabar, X1, weight,Y)
    g1<-t(weight*X1)%*%(weight*Y)/n
    
    if(is.null(L2all)){
      B <- fit1$hessian
      atilde <--(g1+ L1all - fit1$gradient)
    }else{
      B <- L2all
      atilde <- -(g1 + L1all - fit1$gradient)-(L2all-fit1$hessian)%*%betabar
    }
    
    lam.max <- max(abs((atilde + (B - diag(diag(B)))%*%betainit)/diag(B))[-1])
    #lam.min <- 0.02*lam.max
    lam.min <- 0.01*lam.max
    lam.seq <- exp(seq(log(lam.min), log(lam.max), length =nlambda )) 
    
    if(method=="CV"){
      folds=cv.folds(n,nfolds)
      re<-NULL
      for (lambda in lam.seq) {
        se<-0
        for (k in 1:nfolds) {
          train<-as.vector(unlist(folds[-k]))
          Y.train<-Y[train]
          X1.train<-X1[train,]
          weight.train<-weight[train]
          
          test<-as.vector(unlist(folds[k]))
          Y.test<-Y[test]
          X1.test<-X1[test,]
          weight.test<-weight[test]
          
          fit<-adap_ls_ho(betainit, 
                          betabar, 
                          X1=X1.train,Y=Y.train,
                          L1all, L2all, 
                          weight = weight.train,
                          lambda)
          
          beta<-fit$Beta_est_adap
          se=se+sum((weight.test*(Y.test-X1.test%*%beta))^2)
        }
        re<-c(re,se)
      }
      lambda<-lam.seq[which.min(re)]
      fit<-adap_ls_ho(betainit, 
                      betabar, 
                      X1=X1,Y=Y,
                      L1all, L2all, 
                      weight = weight,
                      lambda)
      beta<-fit$Beta_est_adap
      msg<-fit$msg
    }
    
    if(method=="BIC"){
      re<-NULL
      for(lambda in lam.seq){
        fit<-adap_ls_ho(betainit, 
                        betabar, 
                        X1=X1,Y=Y,
                        L1all, L2all, 
                        weight = weight,
                        lambda)
        beta<-fit$Beta_est_adap
        re<-c(re, mean((weight*(Y-X1%*%beta))^2)+log(n)*sum(I(beta!=0)/n))
      }
      lambda<-lam.seq[which.min(re)]
      fit<-adap_ls_ho(betainit, 
                      betabar, 
                      X1=X1,Y=Y,
                      L1all, L2all, 
                      weight = weight,
                      lambda)
      beta<-fit$Beta_est_adap
      msg<-fit$msg
    }
    if(msg==0){
      cat("Successful converge!", "\n")
    }else{
      stop("Lambda is too small!", "\n")
    } 
    
    return(list(coeff=beta, msg=msg))  
  } 
  
  #coordinate descent algorithm for one-shot
  adap_ls_ho<-function(betainit, 
                       betabar,
                       X1, 
                       Y, 
                       L1all, 
                       L2all=NULL,
                       weight,
                       lambda){
    n<-length(Y)
    fit1<-GH_ls_ho( betabar, X1, weight,Y)
    g1<-t(weight*X1)%*%(weight*Y)/n
    
    if(is.null(L2all)){
      B <- fit1$hessian
      atilde <--(g1+ L1all - fit1$gradient)
    }else{
      B <- L2all
      atilde <- -(g1 + L1all - fit1$gradient)-(L2all-fit1$hessian)%*%betabar
    }
    fitcd <- coordi.c(atilde, B, betainit, lambda = lambda, iter_max=4000)
    beta.new <- fitcd$betainit
    
    return(list(Beta_est_adap=beta.new, msg=fitcd$message))
  }
  
  #coordinate descent algorithm
  coordi.c<-function(atilde, B, betainit, lambda, iter_max=100){
    msg<-1
    fit<-coordi_cho(atilde, B, betainit, lambda,iter_max)
    if (fit$iter<= iter_max) {
      #cat("At INNER loop", loop, "Successful converge", "\n")
      msg <- 0
    } 
    return(list(betainit=fit$betainit, message=msg))
  } 
  
  # download derivatives of other sites from the cloud
  # calculate 2nd order approx of the total logL
  logL_all_D1 <- rep(0, px)
  logL_all_D2 <- matrix(0, px, px)
  N <- 0
  for(site_i in control$sites){
    derivatives_i <- pdaGet(paste0(site_i,'_OMderive'),config)
    logL_all_D1 <- logL_all_D1 + derivatives_i$logL_D1*derivatives_i$site_size
    logL_all_D2 <- logL_all_D2 + derivatives_i$logL_D2*derivatives_i$site_size
    N <- N + derivatives_i$site_size
  }
  
  logL_all_D1<-logL_all_D1/N
  logL_all_D2<-logL_all_D2/N
  
  # initial beta
  init_i <- pdaGet(paste0(control$lead_site,'_OMinitialize'),config)
  bbar <- init_i$bhat_i 
  #bbar <- control$beta_init  # derivatives_i$b_meta
  
  fit<-cv.adap_ls_ho( bbar, 
                      bbar,  
                      X1=X,  
                      Y=status, 
                      logL_all_D1, 
                      logL_all_D2,
                      nlambda=20,
                      weight = ipw_weight,
                      nfolds=5,
                      method="BIC")
  
  # Htilde = sol$hessian, 
  surr <- list(btilde = fit$coeff, 
               #setilde=sqrt(diag(solve(sol$hessian))/N), 
               site=config$site_id, 
               site_size=nrow(ipdata))
  
  return(surr)
}





#' @useDynLib pda
#' @title DisC2o AIPW estimate of the ATE at each site
#' 
#' @usage DisC2o.AIPWestimate(ipdata,control,config)
#' @param ipdata local data in data frame
#' @param control pda control
#' @param config pda cloud configuration
#'
#' @return  list(btilde=btilde,  Vtilde=Vtilde)
#' @keywords internal
DisC2o.AIPWestimate <- function(ipdata,control,config) {
  
  status <- ipdata$status
  X <- as.matrix(ipdata[,-c(1,2)])
  treat<-ipdata$treatment
  
  expit = function(x){1/(1+exp(-x))}
  
  AIPW_site<-function(thetahat,
                      betahat,
                      X,
                      Treat,
                      Y){
      mean0<-X%*%thetahat
      weight0<-Treat/expit(mean0)
      mean1<-X%*%betahat
      tau<-mean(mean1+weight0*(Y-mean1))
      V<-mean( (mean1-tau)^2+(weight0^2)*((Y-mean1)^2))
    return(list(tau=tau, V=V)) 
  }
  
  init_i <- pdaGet(paste0(control$lead_site,'_PSestimate'),config)
  thetahat<-init_i$thetatilde
  #thetahat<-control$thetahat
  
  init_i <- pdaGet(paste0(control$lead_site,'_OMestimate'),config)
  betahat<-init_i$btilde
  #betahat<-control$betahat
  
  fit<-AIPW_site(thetahat,betahat,X,treat,status) 
  # message("all surrogate estimates synthesized, no need to broadcast! ")
  return(list(AIPW=fit$tau, V_AIPW=fit$V, site_size=nrow(ipdata)))
}





#' @useDynLib pda
#' @title DisC2o AIPW estimate of the ATE, synthesizing all sites
#' 
#' @usage DisC2o.synthesize(ipdata,control,config)
#' @param ipdata local data in data frame
#' @param control pda control
#' @param config pda cloud configuration
#'
#' @return  list(btilde=btilde,  Vtilde=Vtilde)
#' @keywords internal
DisC2o.synthesize <- function(ipdata,control,config) {

  K <- length(control$sites)
  tau<-0
  V_tau<-0
  N<-0
  for(site_i in control$sites){
    surr_i <- pdaGet(paste0(site_i,'_AIPWestimate'),config)
    nsite<-surr_i$site_size
    N<-N+nsite
    tau<-tau+nsite*surr_i$AIPW
    V_tau<-V_tau+nsite*surr_i$V_AIPW
  }
  tau<-tau/N
  V_tau<-V_tau/N
  
  message("all surrogate estimates synthesized, no need to broadcast! ")
  return(list(AIPW=tau, V_AIPW=V_tau))
}
