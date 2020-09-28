# https://style.tidyverse.org/functions.html#naming

ODAC.steps<-c('initialize','derive','derive_UWZ','estimate','synthesize')
ODAC.family<-'cox'

#' @useDynLib pda
#' @title ODAC initialize
#' 
#' @usage ODAC.initialize(ipdata, control, config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#' @NoRd
#' 
#' @references  Rui Duan, et al. "Learning from local to global: An efficient distributed algorithm for modeling time-to-event data". 
#'               Journal of the American Medical Informatics Association, 2020, https://doi.org/10.1093/jamia/ocaa044
#' @return  list(T_i = T_i, bhat_i = fit_i$coef, Vhat_i = summary(fit_i)$coef[,2]^2, site=control$mysite, site_size= nrow(ipdata))
ODAC.initialize <- function(ipdata,control,config){
    T_i <- sort(unique(ipdata$time[ipdata$status==TRUE]))
    fit_i <- survival::coxph(survival::Surv(time, status) ~ ., data=ipdata)
    
    init <- list(T_i = T_i,
                 bhat_i = fit_i$coef,
                 Vhat_i = summary(fit_i)$coef[,2]^2,   # not as glm, coxph summary can keep NA's! but vcov fills 0's!  
                 site = config$site_id,
                 site_size = nrow(ipdata))
  return(init)
}



#' @useDynLib pda
#' @title Generate pda derivatives
#' 
#' @usage ODAC.derive(ipdata, control, config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#' @NoRd
#' @import Rcpp, RcppArmadillo
#' 
#' @return  list(T_all=T_all, b_meta=b_meta, site=control$mysite, site_size = nrow(ipdata), U=U, W=W, Z=Z, logL_D1=logL_D1, logL_D2=logL_D2)
ODAC.derive <- function(ipdata,control,config) {
    px <- ncol(ipdata) - 2
    # decide if doing ODAC derivatives 1st substep (calculate summary stats U, W, Z) 
    # or 2nd substep (calculate derivatives logL_D1, logL_D2)
    # collect event time pts and meta est from the cloud
    T_all <- c()
    bhat_wt_sum <- rep(0, px)
    wt_sum <- rep(0, px)     # cov matrix?
    for(site_i in control$sites){
      init_i <- pdaGet(paste0(site_i,'_initialize'),config)
      T_all <- c(T_all, init_i$T_i)
      bhat_wt_sum <- bhat_wt_sum + init_i$bhat_i / init_i$Vhat_i
      wt_sum <- wt_sum + 1 / init_i$Vhat_i  # cov matrix?
    }

    T_all <- sort(unique(T_all))
    nt <- length(T_all)
    b_meta <- bhat_wt_sum / wt_sum
    bbar <- b_meta
    # add fake data points to help calculate the summary stats in risk sets ar each time pts
    t_max <- max(ipdata$time)+1
    #generate dataframe in format expected by ODAC
    pfdata <- cbind(T_all, 0, matrix(0, nt, px))
    pfdata <- rbind(ipdata, pfdata, use.names=FALSE)
    pfdata <- pfdata[, 'interval':=cut(pfdata$time, breaks = c(T_all, t_max), labels = 1:nt, right=FALSE)][order(pfdata$interval),]
    pfdata$interval[is.na(pfdata$interval)]<-nt
    X <- as.matrix(pfdata[, control$variables, with=F])
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
    U <- UWZ[,1] - c(nt:1)
    W <- UWZ[,2:(px+1)]
    Z <- array(UWZ[,-c(1:(px+1))], c(nt,px,px))
      
    # summary_stat
    derivatives <- list(T_all=T_all, b_meta=b_meta, site=config$site_id, site_size=nrow(ipdata), U=U, W=W, Z=Z)

  # broadcast to the cloud?
  return(derivatives)
}


#' @useDynLib pda
#' @title Generate pda UWZ derivatives
#' 
#' @usage ODAC.derive_UWZ(ipdata, control, config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#' @NoRd
#' 
#' @details Calculate and broadcast 1st and 2nd order derivative at initial bbar
#'        for ODAC, this requires 2 substeps: 1st calculate summary stats (U, W, Z), 
#'        2nd calculate derivatives (logL_D1, logL_D2)
#'
#' @import Rcpp, RcppArmadillo
#' @return  list(T_all=T_all, b_meta=b_meta, site=control$mysite, site_size = nrow(ipdata), U=U, W=W, Z=Z, logL_D1=logL_D1, logL_D2=logL_D2)
ODAC.derive_UWZ <- function(ipdata,control,config){
    px <- ncol(ipdata) - 2
    # decide if doing ODAC derivatives 1st substep (calculate summary stats U, W, Z) 
    # or 2nd substep (calculate derivatives logL_D1, logL_D2)
      # read and add up (U W Z) from all sites from the cloud
      for(site_i in control$sites){
        sumstat_i <- pdaGet(paste0(site_i,'_derive'),config)
        
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
      d <- c(table(c(ipdata[ipdata$status==T,time], T_all)) - 1)
      
      # 1st and 2nd derivatives
      # X <- as.matrix(ipdata[status==TRUE, control$risk_factor, with=F])
      X <- as.matrix(ipdata[ipdata$status==TRUE, control$variables, with=F])
      
      logL_D1 <- apply(X, 2, sum) - apply(d * W / U, 2, sum, na.rm=T)
      W2 <- array(NA, c(dim(W), px))
      for(ii in 1:px) W2[,,ii] <- W[,ii] * W
      logL_D2 <- apply(d * (W2 - U*Z) / U^2, c(2, 3), sum, na.rm=T)  
      derivatives <- list(T_all=T_all, b_meta=sumstat_i$b_meta, U=U, W=W, Z=Z, 
                          site=config$site_id, site_size = nrow(ipdata),
                          logL_D1=logL_D1, logL_D2=logL_D2)
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
#' @NoRd
#' 
#' @details step-4: construct and solve surrogate logL at the master/lead site
#' @import Rcpp, RcppArmadillo
#' @return  list(btilde = sol$par, Htilde = sol$hessian, site=control$mysite, site_size=nrow(ipdata))
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
      derivatives_i <- pdaGet(paste0(site_i,'_derive_UWZ'),config)
      logL_all_D1 <- logL_all_D1 + derivatives_i$logL_D1
      logL_all_D2 <- logL_all_D2 + derivatives_i$logL_D2
      N <- N + derivatives_i$site_size
    }
    
    # initial beta
    bbar <- derivatives_i$b_meta
    
    # logL at local site
    if(hasTies){
      # rcpp function is negative logL...
      logL_local <- function(beta) 0-rcpp_coxph_logL_efron(beta, time = time, event = status, z = X) # / n
      logL_local_D1 <- function(beta) 0-rcpp_coxph_logL_gradient_efron(beta, time = time, event = status, z = X) # / n
      logL_local_D2 <- function(beta) 0-matrix(rcpp_coxph_logL_hessian(beta, time = time, event = status, z = X), px, px) # / n
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
                      control = list(maxit=control$optim_maxit))

    surr <- list(btilde = sol$par, Htilde = sol$hessian, site=config$site_id, site_size=nrow(ipdata))

    return(surr)
}



#' @useDynLib pda
#' @title PDA synthesize surrogate estimates from all sites, optional
#' 
#' @usage ODAC.synthesize(ipdata, control, config)
#' @param ipdata local data in data frame
#' @param control pda control
#' @param config cloud config
#' @NoRd
#' 
#' @details Optional step-4: synthesize all the surrogate est btilde_i from each site, if step-3 from all sites is broadcasted
#' @import Rcpp, RcppArmadillo
#' @return  list(btilde=btilde,  Vtilde=Vtilde)
ODAC.synthesize <- function(ipdata,control,config) {
  
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
  
  cat("all surrogate estimates synthesized, no need to broadcast! \n")
  return(list(btilde=btilde, 
              Vtilde=Vtilde))
}
