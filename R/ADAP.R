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
#********* remove "synthesize" ***************#
ADAP.steps <- c('initialize','derive','estimate')
ADAP.family <- 'lasso'

#' @useDynLib pda
#' @title ADAP initialize
#' 
#' @usage ADAP.initialize(ipdata,control,config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#' 
#' @return init
#' @keywords internal
ADAP.initialize <- function(ipdata,control,config){
  #****** ipdata columns: outcome, intercept, covariates
  fit0 <- glmnet::cv.glmnet(as.matrix(ipdata[,-c(1:2)]), ipdata$status, family = "binomial")
  beta_i <- as.numeric(coef(fit0, s = "lambda.min"))
  
  #ipdata include intercept, so use 0 + . in glm
  #fit_i <- glm(status ~ 0+., data=ipdata,family = "binomial"(link = "logit"))
  init <- list(site = config$site_id,
               site_size = nrow(ipdata),
               bhat_i = beta_i, #fit_i$coef,
               #**************** changed *********************#
               #use rep to avoid problem in pdasync in line524
               Vhat_i = rep(nrow(ipdata),length(beta_i)))
  return(init)
}

#' @useDynLib pda
#' @title ADAP derivatives
#' 
#' @usage ADAP.derive(ipdata,control,config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#'
#' @return  list(site=config$site_id, site_size = nrow(ipdata), logL_D1=logL_D1, logL_D2=logL_D2)
#' @keywords internal
ADAP.derive <- function(ipdata,control,config){
  # data sanity check ...
  px <- ncol(ipdata) - 1  # X includes intercept
  # get b_meta as initial bbar
  bhat <- rep(0, px)
  vbhat <- rep(0, px)  # sample size to be used as weights in ADAP  
  #******** control$sites contains all sites
  for(site_i in control$sites){
    init_i <- pdaGet(paste0(site_i,'_initialize'),config)
    bhat = rbind(bhat, init_i$bhat_i)
    vbhat = rbind(vbhat, init_i$Vhat_i)
  }
  bhat = bhat[-1,]
  vbhat = vbhat[-1,]
  
  # #estimate from meta-analysis
  # betameta = apply(bhat/vbhat,2,function(x){sum(x, na.rm = T)})/apply(1/vbhat,2,function(x){sum(x, na.rm = T)})
  # vmeta = 1/apply(1/vbhat,2,function(x){sum(x, na.rm = T)})
  
  # initialize with the weighted average 
  betameta <- apply(diag(vbhat[,1]/sum(vbhat[,1]))%*%bhat, 2, sum, na.rm = TRUE)
  bbar <- betameta 
  
  # 1st and 2nd derivatives
  status <- ipdata$status
  X <- as.matrix(ipdata[,-1])
  
  expit = function(x){1/(1+exp(-x))}
  
  #******** Rui: logL, ADAP: -logL as objective function
  #first order gradient
  Lgradient = function(beta,X,Y){
    design = X  # cbind(1,X)
    -t(Y-expit(design%*%t(t(beta))))%*%design/length(Y)
  }
  
  #second-order gradient
  Lgradient2 = function(beta,X){
    design = X # cbind(1,X)
    Z=expit(design%*%beta)
    t(c(Z*(1-Z))*design)%*%design/nrow(X)
  }
  
  logL_D1 <- Lgradient(bbar,X,ipdata$status)
  logL_D2 <- Lgradient2(bbar,X)
  
  derivatives <- list(
    site=config$site_id, 
    site_size = nrow(ipdata),
    logL_D1=logL_D1,
    logL_D2=logL_D2)  # , bbar=bbar, bhat=bhat
  
  return(derivatives)
}


#' @useDynLib pda
#' @title ADAP surrogate estimation
#' 
#' @usage ADAP.estimate(ipdata,control,config)
#' @param ipdata local data in data frame
#' @param control PDA control
#' @param config cloud configuration
#' 
#' @details step-3: construct and solve surrogate objective function at the master/lead site
#' @return  list(btilde = sol$par, Htilde = sol$hessian, site=control$mysite, site_size=nrow(ipdata))
#' @keywords internal
ADAP.estimate <- function(ipdata,control,config) {
  # data sanity check ...
  status <- ipdata$status
  #******** X includes 1
  X <- as.matrix(ipdata[,-1])
  px <- ncol(X)
  n <- nrow(X)
  Y <- ipdata$status
  
  ######################################################
  # #likelihood function for logistic regression, the input X is a n*d matrix where
  # #each patient has d covariates stored in each row.
  # Lik = function(beta,X,Y){
  #   design = X # cbind(1,X)
  #   sum(Y*(design%*%t(t(beta)))-log(1+exp(design%*%t(t(beta)))))/length(Y)
  # }
  
  expit = function(x){1/(1+exp(-x))}
  
  #loss function: negative log-likelihood function for logistic regression, the input X is a n*d matrix
  NLogLik <- function(beta, X, Y){
    design = X #cbind(1, X)
    -sum(Y*(design%*%t(t(beta))) - log(1 + exp(design%*%t(t(beta)))))/length(Y)
  }
  
  #******** Rui: logL, ADAP: -logL as objective function
  #first order gradient
  Lgradient = function(beta,X,Y){
    design = X  # cbind(1,X)
    -t(Y-expit(design%*%t(t(beta))))%*%design/length(Y)
  }
  
  #second-order gradient
  Lgradient2 = function(beta,X){
    design = X # cbind(1,X)
    Z=expit(design%*%beta)
    t(c(Z*(1-Z))*design)%*%design/nrow(X)
  }
  
  #soft-thresholding operator
  soft <- function(x,thres){
    if (x > thres){
      res <- x - thres
    } else if (x < -thres){
      res <- x + thres
    } else {
      res <- 0
    }
    res
  }
  
  #coordinate descent algorithm
  coordi <- function(atilde, B, betainit, lambda){
    tol <- 1e-4
    loop <- 1
    Loop <- 100
    msg <- NA
    p <- length(atilde)
    # if (any(diag(B)) == 0) {
    #   cat("Error: exist zero variable", "\n")
    #   msg <- "Error"
    #   break
    # }
    while (loop <= Loop) {
      beta.old <- betainit
      for (j in 1:p){
        a <- diag(B)[j]
        ej <- diag(p)[j,]
        b <- atilde[j] + ej%*%B%*%betainit - a*betainit[j]
        if (j == 1) {
          #no penalization on intercept
          betainit[1] <- -b/a 
        } else {
          betainit[j] <- -soft(b/a,lambda/a)
        }
      }
      # dif <- sum((betainit-beta.old)^2)
      # use the change in loss function to measure 
      dif <- (atilde)%*%(betainit-beta.old) + (betainit)%*%B%*%(betainit)/2 - 
        (beta.old)%*%B%*%(beta.old)/2
      if (abs(dif) < tol) {
        cat("At INNER loop", loop, "Successful converge", "\n")
        msg <- "Successful convergence"
        break
      }
      if (loop == Loop) {
        cat("INNER loop maximum iteration reached", "\n")
      } 
      loop <- loop + 1
    }
    list(betainit=betainit, message=msg)
  }
  
  
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
  
  # initial beta use the averge estimator
  bbar <- c(control$beta_init) # derivatives_i$bbar   # 
  # bhat <- derivatives_i$bhat
  # use the 1st site to expand the objective function
  # betatilde <- bhat[1,]     
  betatilde <- c(pdaGet(paste0(control$lead_site,'_initialize'),config)$bhat_i )
  # print(length(bbar))
  # print(length(betatilde))
  
  
  # calculate the gradient at betatilde
  Ltilde <- Lgradient(betatilde, X, Y)
  L2tilde <- as.vector(Lgradient2(betatilde, X))
  
  # Calculate the global first order gradient L
  L_all <- logL_all_D1/N
  # Calculate the global second order gradient L2
  L2_all <- as.vector(logL_all_D2/N)  #  why as.vector and then matrix(L2_all, ncol = p, nrow = p) ? could be slow...
  print(dim(t(L_all)))
  print(dim(logL_all_D2))
  
  Beta_est <- rep(NA, px)
  Xlocal <- X
  Ylocal <- Y
  p <- px
  betabar <- bbar
  
  ###############
  ###  ODAL2  ###
  ###############
  
  #use betatilde to get quadratic approximation
  B <- Lgradient2(betatilde, Xlocal) + matrix(L2_all, ncol = p, nrow = p) - Lgradient2(betabar, Xlocal) 
  tt = Lgradient2(betatilde, Xlocal)
  print(dim(t(betatilde)%*%Lgradient2(betatilde, Xlocal)))
  atilde <- Lgradient(betatilde, Xlocal, Ylocal) - t(betatilde)%*%Lgradient2(betatilde, Xlocal) + L_all - 
            Lgradient(betabar, Xlocal, Ylocal) - t(betabar)%*%(matrix(L2_all, ncol = p, nrow = p) - Lgradient2(betabar, Xlocal))
  
  lam.max <- max(abs(atilde + t((B - diag(diag(B)))%*%betabar)))
  lam.min <- ifelse(n < p, 0.02, 1e-04)*lam.max
  lam.seq <- exp(seq(log(lam.min), log(lam.max), length = 100))
  #lam.seq <- rev(cv.glmnet(Xlocal, Ylocal, family = "binomial", nfolds = 5)$lambda)
  
  ##########################
  ###  5CV  ###
  ##########################
  
  nfold <- 5
  norder <- NULL
  if (is.null(norder)) norder <- sample(seq_len(n),n)
  
  lam_path <- matrix(ncol = nfold, nrow = length(lam.seq), NA)
  
  ndel <- round(n/nfold)
  for (f in seq_len(nfold)){
    if (f != nfold) {
      iddel <- norder[(1 + ndel * (f - 1)):(ndel * f)]
    } else {
      iddel <- norder[(1 + ndel * (f - 1)):n]
    }
    ndel <- length(iddel)
    nf <- n - ndel
    idkeep <- (seq_len(n))[-iddel]
    
    Xf <- Xlocal[-iddel, ]
    Xfdel <- Xlocal[iddel, ]
    Yf <- as.matrix(Ylocal[-iddel])
    Yfdel <- as.matrix(Ylocal[iddel])
    
    # adjust for Cross Validation procedure
    # the global function excludes the validation set 
    L2_allf <- (L2_all*N - as.vector(Lgradient2(betabar, Xfdel))*ndel)/(N - ndel)
    Bf <- Lgradient2(betatilde, Xf) + matrix(L2_allf, ncol = p, nrow = p) - 
      Lgradient2(betabar, Xf) 
    atilde_f <- Lgradient(betatilde, Xf, Yf) - t(betatilde)%*%Lgradient2(betatilde, Xf) + 
      (L_all*N - Lgradient(betabar, Xfdel, Yfdel)*ndel)/(N - ndel) - 
      Lgradient(betabar, Xf, Yf) - 
      t(betabar)%*%(matrix(L2_allf, ncol = p, nrow = p) - Lgradient2(betabar, Xf))
    betainit <- rep(0, p)
    
    for (la in 1:length(lam.seq)){
      tol <- 1e-5
      out.loop <- 1
      out.Loop <- 100
      msg <- 0
      
      while (out.loop <= out.Loop) {
        fitcd <- coordi(atilde_f, Bf, betainit, lambda = rev(lam.seq)[la])
        beta.new <- fitcd$betainit
        
        if (is.na(fitcd$message)){
          msg <- 1
          cat("At OUTER loop", out.loop, "Lambda is too small", "\n")
          break
        }
        #dif <- sum((beta.new-betainit)^2)
        dif <- (atilde_f)%*%(beta.new - betainit) + (beta.new)%*%Bf%*%(beta.new)/2 - 
          (betainit)%*%Bf%*%(betainit)/2
        
        Bf <- Lgradient2(beta.new, Xf) + matrix(L2_allf, ncol = p, nrow = p) - 
          Lgradient2(betabar, Xf) 
        atilde_f <- Lgradient(beta.new, Xf, Yf) - t(beta.new)%*%Lgradient2(beta.new, Xf) + 
          (L_all*N - Lgradient(betabar, Xfdel, Yfdel)*ndel)/(N - ndel) - 
          Lgradient(betabar, Xf, Yf) - 
          t(betabar)%*%(matrix(L2_allf, ncol = p, nrow = p) - Lgradient2(betabar, Xf))
        
        betainit <- beta.new
        if (abs(dif) < tol) {
          cat("At OUTER loop", out.loop, "Successful converge", "\n")
          break
        }
        if (out.loop == out.Loop) {
          cat("OUTER loop maximum iteration reached", "\n")
          break
        } 
        out.loop <- out.loop + 1
      }
      if (msg==1) break
      if (out.loop < out.Loop) lam_path[la,f] <- NLogLik(beta.new, Xfdel, Yfdel)
    }
  }
  index <- order(colSums(lam_path))
  crerr <- rowSums(lam_path[, index])/length(index) * nfold
  lam.est <- rev(lam.seq)[which.min(crerr)]
  
  #fit the final estimator
  betainit <- betatilde
  tol <- 1e-5
  out.loop <- 1
  out.Loop <- 100
  msg <- 0
  while (out.loop <= out.Loop) {
    fitcd <- coordi(atilde, B, betainit, lambda = lam.est)
    if (is.na(fitcd$message)){
      cat("The local lambda is too small", "\n")
      msg <- 1
      break
    }
    beta.new <- fitcd$betainit
    
    #dif <- sum((beta.new-betainit)^2)
    dif <- (atilde)%*%(beta.new - betainit) + (beta.new)%*%B%*%(beta.new)/2 - 
      (betainit)%*%B%*%(betainit)/2
    B <- Lgradient2(beta.new, Xlocal) + matrix(L2_all, ncol = p, nrow = p) - 
      Lgradient2(betabar, Xlocal) 
    atilde <- Lgradient(beta.new, Xlocal, Ylocal) - t(beta.new)%*%Lgradient2(beta.new, Xlocal) + 
      L_all - Lgradient(betabar, Xlocal, Ylocal) -
      t(betabar)%*%(matrix(L2_all, ncol = p, nrow = p) - Lgradient2(betabar, Xlocal))
    
    betainit <- beta.new
    if (abs(dif) < tol) {
      cat("At OUTER loop", out.loop, "Successful converge", "\n")
      break
    }
    if (out.loop == out.Loop) {
      cat("OUTER loop maximum iteration reached", "\n")
      break
    }
    out.loop <- out.loop + 1
  }
  if ((msg != 1)&(out.loop < out.Loop)) Beta_est <- beta.new
  
  
  # #second-order surogate likelihood, suppose the local data are stored in Xlocal, Ylocal
  # Y = ipdata$status
  # n1 = length(Y)
  # logL_tilde = function(beta){
  #   - (Lik(beta,X, Y) + (logL_all_D1/N - Lgradient(bbar, X, Y))%*%beta+
  #        t(beta-bbar)%*%(logL_all_D2/N - Lgradient2(bbar, X))%*%(beta-bbar) / 2)
  # }
  # 
  # # optimize the surrogate logL
  # sol <- optim(par = bbar,
  #              fn = logL_tilde,
  #              # gr = logL_tilde_D1,
  #              hessian = TRUE,
  #              control = list(maxit=control$optim_maxit))
  
  surr <- list(btilde = Beta_est, 
               Htilde = NA, # sol$hessian, 
               site=config$site_id, 
               site_size=nrow(ipdata))
  ######################################################
  
  return(surr)
}
