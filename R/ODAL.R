# https://style.tidyverse.org/functions.html#naming

ODAL.steps<-c('initialize','derive','estimate','synthesize')
ODAL.family<-'binomial'
# require(survival)
# require(data.table)
# Rcpp::sourceCpp('src/rcpp_coxph.cpp')
## broadcast: upload/download shared info to/from the cloud folder

#' @useDynLib pda
#' @title PDA initialize
#' 
#' @usage pda_initialize <- function(ipdata, broadcast=TRUE, control=control)
#' @author Chongliang Luo, Steven Vitale
#' @param ipdata local data in data frame
#' @param control pda control
#' @return  list(T_i = T_i, bhat_i = fit_i$coef, Vhat_i = summary(fit_i)$coef[,2]^2, site=control$mysite, site_size= nrow(ipdata))
ODAL.initialize <- function(ipdata,control,config){
  # data sanity check ...
  
  # if(!any(names(ipdata)[1:2] == c('time', 'status')))
  #   error('ipdata columns should be (time, status, covariates)')
  # if(!any(is.numeric(ipdata)))
  #   error('ipdata need to be numeric, please create dummy variables if necessary')
  
    fit_i <- glm(status ~ 0+., data=ipdata,family = "binomial"(link = "logit"))
    init <- list(site = config$site_id,
                 site_size = nrow(ipdata),
                 bhat_i = fit_i$coef,
                 Vhat_i = summary(fit_i)$coef[,2]^2)
  return(init)
}


#' @useDynLib pda
#' @title PDA derive
#' 
#' @usage ODAL.derive(bbar, ipdata, broadcast=TRUE, derivatives_ODAC_substep='first', control=control)
#' @author Chongliang Luo, Steven Vitale
#' 
#' @param bbar  initial estimate
#' @param ipdata local data in data frame
#' @param broadcast Logical, broadcast to the cloud? 
#' @param derivatives_ODAC_substep character, only for Cox regression, 'first' / 'second' indicate which substep of ODAC   
#' @param control PDA control
#' 
#' @details ONly for ODAC: step-2: calculate and broadcast 1st and 2nd order derivative at initial bbar
#'        for ODAC, this requires 2 substeps: 1st calculate summary stats (U, W, Z), 
#'        2nd calculate derivatives (logL_D1, logL_D2)
#'
#' @return  list(T_all=T_all, b_meta=b_meta, site=control$mysite, site_size = nrow(ipdata), U=U, W=W, Z=Z, logL_D1=logL_D1, logL_D2=logL_D2)
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
      # b_init=bbar,
      logL_D1=logL_D1,
      logL_D2=logL_D2)
  
  # broadcast to the cloud?
  return(derivatives)
}


#' @useDynLib pda
#' @title PDA surrogate estimation
#' 
#' @usage pda_surrogate_est(bbar, ipdata, broadcast=TRUE, control=control)
#' @author Chongliang Luo, Steven Vitale
#' 
#' @param bbar  initial estimate
#' @param ipdata local data in data frame
#' @param broadcast Logical, broadcast to the cloud? 
#' @param control PDA control
#' 
#' @details step-3: construct and solve surrogate logL at the master/lead site
#'
#' @return  list(btilde = sol$par, Htilde = sol$hessian, site=control$mysite, site_size=nrow(ipdata))
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
#' @usage pda_synthesize(control=control)
#' @author Chongliang Luo, Steven Vitale
#' 
#' @param control PDA control
#' 
#' @details Optional step-4: synthesize all the surrogate est btilde_i from each site, if step-3 from all sites is broadcasted
#'
#' @return  list(btilde=btilde,  Vtilde=Vtilde)
ODAL.synthesize <- function(ipdata,control,config) {
  
  px <- length(control$variable)
  K <- length(control$sites)
  # if (control$model == "ODAL" | control$model == "ODALR"){
    # btilde_wt_sum <- rep(0, px+1)
    # wt_sum <- rep(0, px+1)     # cov matrix?
  # }else{
    btilde_wt_sum <- rep(0, px)
    wt_sum <- rep(0, px)     # cov matrix?
  # }
  
  for(site_i in control$sites){
    surr_i <- pdaGet(paste0(site_i,'_estimate'),config)
    btilde_wt_sum <- btilde_wt_sum + surr_i$Htilde %*% surr_i$btilde
    wt_sum <- wt_sum + surr_i$Htilde
  }
  
  # inv-Var weighted average est, and final Var = average Var-tilde
  btilde <- solve(wt_sum, btilde_wt_sum)
  Vtilde <- solve(wt_sum) * K
  
  cat("all surrogate estimates surr_i$Htildesynthesized, no need to broadcast! \n")
  return(list(btilde=btilde, 
              Vtilde=Vtilde))
}
