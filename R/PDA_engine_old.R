

# https://style.tidyverse.org/functions.html#naming

# require(survival)
# require(data.table)
# Rcpp::sourceCpp('src/rcpp_coxph.cpp')

## broadcast: upload/download shared info to/from the cloud folder
## write to csv or excel files for better observation on the cloud?
pda_broadcast <- function(obj,                    # list: object to be broadcasted
                          obj_type=c('initialize', 'summary_stat', 'derivatives', 'surrogate_est'),
                          upload=TRUE,
                          site_i,                 # name of site on the cloud, only for upload==FALSE (i.e. download from cloud)
                          control){
  
  ff <- paste0(control$cloud, '/', control$outcome, '_', control$model, '_', 
               ifelse(upload==TRUE, control$local_site, site_i), '_', obj_type, '.rds')
  if(upload==TRUE)
    saveRDS(obj, file=ff)
  else
    return(readRDS(ff))
}


## step-1: initialize, such as (bhat, Vhat) from each site
## output: list(T_i = T_i, bhat_i = fit_i$coef, Vhat_i = summary(fit_i)$coef[,2]^2, site=control$local_site, site_size= nrow(mydata))
pda_initialize <- function(mydata,          # dataframe including: time, status and X (numeric)
                           broadcast=TRUE,
                           control=pda_control){
  # data sanity check ...
  
  # if(!any(names(mydata)[1:2] == c('time', 'status')))
  #   error('mydata columns should be (time, status, covariates)')
  # if(!any(is.numeric(mydata)))
  #   error('mydata need to be numeric, please create dummy variables if necessary')
  
  if(control$model=='ODAC'){
    T_i <- sort(unique(mydata$time[mydata$status==TRUE]))
    fit_i <- coxph(Surv(time, status) ~ ., data=mydata)
    
    init <- list(T_i = T_i,
                 bhat_i = fit_i$coef,
                 Vhat_i = summary(fit_i)$coef[,2]^2,   # cov matrix? vcov(fit_i)
                 site = control$local_site,
                 site_size = nrow(mydata))
  }

   
  if(control$model=='ODACH'){
    # any diagnosis of heterogeneous baseline hazard?
    
    fit_i <- coxph(Surv(time, status) ~ ., data=mydata)
    
    init <- list(bhat_i = fit_i$coef,
                 Vhat_i = summary(fit_i)$coef[,2]^2,   # cov matrix? vcov(fit_i)
                 site = control$local_site,
                 site_size = nrow(mydata))
  }
  
  
  if(control$model=='ODAL'){
    
  }
  
  
  if(control$model=='ODALH'){
    
  }
  
  
  if(broadcast){
    # if(not exist ...)
    pda_broadcast(init, 'initialize', control=control)
  }else{
    cat('local initialization not broadcasted. Please review the output and 
        use pda_broadcast() to manually upload to the cloud!')
  }
  
  return(init)
}



## step-2: calculate and broadcast 1st and 2nd order derivative at initial bbar
##        for ODAC, this requires 2 substeps: 1st calculate summary stats (U, W, Z), 
##        2nd calculate derivatives (logL_D1, logL_D2)
## output: list(T_all=T_all, b_meta=b_meta, site=control$local_site, site_size = nrow(mydata),  
##              U=U, W=W, Z=Z,
##              logL_D1=logL_D1, logL_D2=logL_D2)
pda_derivatives <- function(bbar = NULL,
                             mydata, 
                             broadcast = TRUE,
                             derivatives_ODAC_substep='first',  # 'second'
                             control = pda_control){
  # data sanity check ...
  
  px <- ncol(mydata) - 2
  
  if(control$model=='ODAC'){
    # decide if doing ODAC derivatives 1st substep (calculate summary stats U, W, Z) 
    # or 2nd substep (calculate derivatives logL_D1, logL_D2)
    if(is.null(derivatives_ODAC_substep)){
      if(any(grepl('derivatives_UWZ', list.files(control$cloud))))
        derivatives_ODAC_substep <- 'second'
      else
        derivatives_ODAC_substep <- 'first'
    }
    
    if(derivatives_ODAC_substep == 'first'){
      # collect event time pts and meta est from the cloud
      T_all <- c()
      bhat_wt_sum <- rep(0, px)
      wt_sum <- rep(0, px)     # cov matrix?
      for(site_i in control$all_site){
        init_i <- pda_broadcast(obj_type= 'initialize',
                                upload=FALSE,
                                site_i=site_i, 
                                control=control) 
        T_all <- c(T_all, init_i$T_i)
        bhat_wt_sum <- bhat_wt_sum + init_i$bhat_i / init_i$Vhat_i
        wt_sum <- wt_sum + 1 / init_i$Vhat_i  # cov matrix?
      }
      
      T_all <- sort(unique(T_all))
      nt <- length(T_all)
      b_meta <- bhat_wt_sum / wt_sum
      if(is.null(bbar)) bbar <- b_meta
      
      # add fake data points to help calculate the summary stats in risk sets ar each time pts
      t_max <- max(mydata$time)+1
      tmp <- cbind(T_all, 0, matrix(0, nt, px))
      tmp <- rbind(mydata, tmp, use.names=FALSE)
      tmp <- tmp[, interval:=cut(time, breaks = c(T_all, t_max), labels = 1:nt, right=F)][order(interval),]
      X <- as.matrix(tmp[, control$risk_factor, with=F])
      
      # summary stats: U, W, Z
      eXb <- c(exp(X %*% bbar))
      X2 <- X[,1]*X
      for(ix in 2:ncol(X)) X2 <- cbind(X2, X[,ix]*X)
      UWZ <- eXb * cbind(1, X, X2)
      
      # rcpp_aggregate() is a function written in rcpp for calculating column-wise (reverse) cumsum
      # credit to Dr Wenjie Wang
      UWZ <- rcpp_aggregate(x = UWZ, indices = tmp$interval, cumulative = T, reversely = T)
      
      # since fake X=0, cumulative W and Z will be the same, 
      # but exp(Xb)=1, so need to remove cumulated ones from each time pts
      U <- UWZ[,1] - c(nt:1)
      W <- UWZ[,2:(px+1)]
      Z <- array(UWZ[,-c(1:(px+1))], c(nt,px,px))
      
      # summary_stat
      derivatives <- list(T_all=T_all, b_meta=b_meta, site=control$local_site, site_size=nrow(mydata), U=U, W=W, Z=Z)
    }
    
    if(derivatives_ODAC_substep == 'second'){
      # read and add up (U W Z) from all sites from the cloud
      for(site_i in control$all_site){
        sumstat_i <- pda_broadcast(obj_type= 'derivatives_UWZ',
                                   upload=FALSE,
                                   site_i=site_i, 
                                   control=control)
        if(site_i == control$all_site[1]){
          U <- sumstat_i$U
          W <- sumstat_i$W
          Z <- sumstat_i$Z
        }else{
          U <- U + sumstat_i$U
          W <- W + sumstat_i$W
          Z <- Z + sumstat_i$Z
        }
      }
      
      # number of events in mydata at each event time pts in T_all
      T_all <- sumstat_i$T_all
      d <- c(table(c(mydata[status==T,time], T_all)) - 1)
    
      # 1st and 2nd derivatives
      X <- as.matrix(mydata[status==TRUE, control$risk_factor, with=F])
      logL_D1 <- apply(X, 2, sum) - apply(d * W / U, 2, sum, na.rm=T)
      W2 <- array(NA, c(dim(W), px))
      for(ii in 1:px) W2[,,ii] <- W[,ii] * W
      logL_D2 <- apply(d * (W2 - U*Z) / U^2, c(2, 3), sum, na.rm=T)  
      
      derivatives <- list(T_all=T_all, b_meta=sumstat_i$b_meta, U=U, W=W, Z=Z, 
                          site=control$local_site, site_size = nrow(mydata),
                          logL_D1=logL_D1, logL_D2=logL_D2)
    }
  }
  
  
  if(control$model=='ODACH'){
    # get b_meta as initial bbar
    bhat_wt_sum <- rep(0, px)
    wt_sum <- rep(0, px)     # cov matrix?
    for(site_i in control$all_site){
      init_i <- pda_broadcast(obj_type= 'initialize',
                              upload=FALSE,
                              site_i=site_i, 
                              control=control) 
      bhat_wt_sum <- bhat_wt_sum + init_i$bhat_i / init_i$Vhat_i
      wt_sum <- wt_sum + 1 / init_i$Vhat_i  # cov matrix?
    }
    
    b_meta <- bhat_wt_sum / wt_sum
    if(is.null(bbar)) bbar <- b_meta
           
    # 1st and 2nd derivatives
    time <- mydata$time
    status <- mydata$status
    X <- as.matrix(mydata[,-c(1,2)])
    n <- length(time)
    px <- ncol(X)
    hasTies <- any(duplicated(mydata$time))
    
    if(hasTies){
      # rcpp function is negative logL...
      logL_D1 <- -rcpp_coxph_logL_gradient_efron(beta = bbar, time = time, event = status, z = X) # / n
      logL_D2 <- -matrix(rcpp_coxph_logL_hessian(beta = bbar, time = time, event = status, z = X), px, px) # / n
    } else {
      logL_D1 <- -rcpp_coxph_logL_gradient(beta = bbar, time = time, event = status, z = X) # / n
      logL_D2 <- -matrix(rcpp_coxph_logL_hessian(beta = bbar, time = time, event = status, z = X), px, px) # / n
    }
    
    derivatives <- list(b_meta=sumstat_i$b_meta,  site=control$local_site, site_size = nrow(mydata),
                        logL_D1=logL_D1, logL_D2=logL_D2)

  }
  
  
  if(control$model=='ODAL'){
    
    
    
    derivatives <- list(#T_all=T_all,  U=U, W=W, Z=Z, 
                        b_meta=sumstat_i$b_meta,
                        site=control$local_site, site_size = nrow(mydata),
                        logL_D1=logL_D1, logL_D2=logL_D2)
  }
  
  
  if(control$model=='ODALH'){
    
  }
  
  
  
  # broadcast to the cloud?
  if(broadcast){
    # if(not exist ...)
    obj_type <- 'derivatives'
    if(control$model=='ODAC' & derivatives_ODAC_substep == 'first')    
      obj_type <- 'derivatives_UWZ'

    pda_broadcast(derivatives, obj_type, control=control)
  }else{
    cat('Derivatives not broadcasted. Please review the output and 
          use pda_broadcast() to manually upload to the cloud!')
  }
  
  return(derivatives)
}



## step-3: construct and solve surrogate logL at my site
## output: list(btilde = sol$par, Htilde = sol$hessian, site=control$local_site, site_size=nrow(mydata))
pda_surrogate_est <- function(bbar = NULL,
                              mydata, 
                              broadcast = FALSE,       # broadcasting is optional for this step
                              control = pda_control){
  # data sanity check ...
  
  if(control$model=='ODAC' | control$model=='ODACH'){
    time <- mydata$time
    status <- mydata$status
    X <- as.matrix(mydata[,-c(1,2)])
    n <- length(time)
    px <- ncol(X)
    hasTies <- any(duplicated(mydata$time))
    
    # download derivatives of other sites from the cloud
    # calculate 2nd order approx of the total logL  
    logL_all_D1 <- rep(0, px)
    logL_all_D2 <- matrix(0, px, px)
    N <- 0
    for(site_i in control$all_site){
      derivatives_i <- pda_broadcast(obj_type= 'derivatives',
                                 upload=FALSE,
                                 site_i=site_i, 
                                 control=control)
      logL_all_D1 <- logL_all_D1 + derivatives_i$logL_D1
      logL_all_D2 <- logL_all_D2 + derivatives_i$logL_D2
      N <- N + derivatives_i$site_size
    }
    
    # initial beta
    if(is.null(bbar)) bbar <- derivatives_i$b_meta
    
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

    surr <- list(btilde = sol$par, Htilde = sol$hessian, site=control$local_site, site_size=nrow(mydata))
  }

   
  
  if(control$model=='ODAL'){
    
    
  
    
    # # optimize the surrogate logL 
    # sol <- optim(par = bbar, 
    #              fn = logL_tilde,
    #              # gr = logL_tilde_D1,
    #              hessian = TRUE, 
    #              control = list(maxit=control$optim_maxit))
    # 
    # surr <- list(btilde = sol$par, Htilde = sol$hessian, site=control$local_site, site_size=nrow(mydata))
  }

    
  if(control$model=='ODALH'){
  
  }
  
  # broadcast to the cloud?
  if(broadcast){
    # if(not exist ...)
    pda_broadcast(surr, 'surrogate_est', control=control)
  }else{
    cat('Surrogate estimate not broadcasted. Please review the output and 
          use pda_broadcast() to manually upload to the cloud!')
  }
  
  return(surr)
}




## step-4: synthesize all the surrogate est btilde_i from each site, if step-3 from all sites is broadcasted
pda_synthesize <- function(control = pda_control){
  
  px <- length(control$risk_factor)
  K <- length(control$all_site)
  btilde_wt_sum <- rep(0, px)
  wt_sum <- rep(0, px)     # cov matrix?
  for(site_i in control$all_site){
    surr_i <- pda_broadcast(obj_type = 'surrogate_est',
                               upload=FALSE,
                               site_i=site_i, 
                               control=control)
    btilde_wt_sum <- btilde_wt_sum + surr_i$Htilde %*% surr_i$btilde
    wt_sum <- wt_sum + surr_i$Htilde
  }
  
  # inv-Var weighted average est, and final Var = average Var-tilde
  btilde <- solve(wt_sum, btilde_wt_sum)
  Vtilde <- solve(wt_sum) * K
  
  return(list(btilde=btilde, 
              Vtilde=Vtilde))
}



pda_main <- function(mydata = mydata,
                     step = 1,                       # c('initialize', 'derivatives', 'surrogate_est', 'synthesize'),
                     derivatives_ODAC_substep=NULL,  # c('first', 'second'),  only for control$model=='ODAC' and step==2
                     control = pda_control){
  if(step==1 | step=='initialize'){
    output <- pda_initialize(mydata, control=control)
    print(output$bhat_i)
    print(output$site)
    print(output$site_size)
  } 
  
  # if(step==2 | step=='summary_stat'){
  #   output <- pda_summary_stat(bbar=NULL, mydata, control=control)
  #   print(output$b_meta) 
  # }
  
  if(step==2 | step=='derivatives') 
    output <- pda_derivatives(bbar=NULL, mydata, derivatives_ODAC_substep=derivatives_ODAC_substep, control=control)
    
  
  if(step==3 | step=='surrogate_est'){
    output <- pda_surrogate_est(bbar=NULL, mydata, control=control)
    cat('\n', output$btilde)
    return(output)
  }
  
  if(step==4 | step=='synthesize'){
    output <- pda_synthesize(control=control)
    print(output$btilde)
    return(output)
  }
  
  # return(output)
}


#