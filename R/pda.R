# https://style.tidyverse.org/functions.html#naming

# require(survival)
# require(data.table)
# Rcpp::sourceCpp('src/rcpp_coxph.cpp')

## broadcast: upload/download shared info to/from the cloud folder
## write to csv or excel files for better observation on the cloud?


#' @useDynLib PDA
#' @title Function to communicate files (upload/download summary statistics to/from the cloud)
#' 
#' @usage pda_broadcast(obj, obj_type=c('initialize', 'summary_stat', 'derivatives', 'surrogate_est'), file_name = NULL, upload=TRUE, site_i, control)
#' @author Chongliang Luo, Steven Vitale
#' 
#' @param obj R object to be broadcasted
#' @param obj_type object type to be communicated, can be 'initialize', 'summary_stat', 'derivatives', 'surrogate_est' depends on which PDA step
#' @param file_name file name to be communicated
#' @param upload   Logical, TRUE/FALSE = upload/download summary statistics to/from the cloud
#' @param site_i name of site on the cloud, only for upload==FALSE (i.e. download from cloud)
#' @param control PDA control
#'
#' @return  
#' @export
pda_broadcast <- function(obj,                    
                          obj_type=c('initialize', 'summary_stat', 'derivatives', 'surrogate_est'),
                          file_name = NULL,
                          upload=TRUE,
                          site_i,                 
                          control){
  
  if(is.null(file_name)){
    site = ifelse(upload==TRUE, control$mysite, site_i)
    file_name = paste0(site, '_', obj_type)
  }
  ff <- paste0(tempdir(), '/' , file_name)
  if(upload==TRUE){
    cat('The following summary statistics have been uploaded to the public cloud: \n')
    print(obj)
    pda_put(obj,file_name)
  } else{
    res = pda_get(file_name)
    cat('The following summary statistics have been read from the public cloud: \n')
    print(res)
    return(res)
  }
    
}

#' @useDynLib PDA
#' @title Function to upload object as RDS
#' 
#' @usage pda_put(obj,name)
#' @author Chongliang Luo, Steven Vitale
#' 
#' @param obj R object to upload
#' @param name of file
#'
#' @return  
#' @export
pda_put <- function(obj,name){
    file_name <- paste0(name, '.json')
    # the file to upload
    file_path <- paste0(tempdir(),'/', file_name)
    obj_Json <- jsonlite::toJSON(obj)
    write(obj_Json, file_path)
    # create the url target of the file
    if (Sys.getenv("PDA_URI") != "") {
        url <- file.path(Sys.getenv("PDA_URI"), file_name)
        # webdav uses a PUT request to send a file to Nextcloud
        r<-httr::PUT(url, body = upload_file(file_path), authenticate(Sys.getenv('PDA_SITE'), Sys.getenv('PDA_SECRET'), 'digest'))
        return(r)
    }
}


#' @useDynLib PDA
#' @title Function to download RDS and return as object)
#' 
#' @usage pda_get(name)
#' @author Chongliang Luo, Steven Vitale
#' 
#' @param name of file
#'
#' @return  
#' @export
pda_get <- function(name){
    print("getting file")
    print(name)
    file_name <- paste0(name, '.json')
    # the file to create
    file_path <- paste0(tempdir(),'/', file_name)
    # create the url target of the file
    if (Sys.getenv("PDA_URI") != "") {
        url <- file.path(Sys.getenv("PDA_URI"), file_name)
        # webdav uses a PUT request to send a file to Nextcloud
        res<-httr::GET(url, write_disk(file_path, overwrite = TRUE), authenticate(Sys.getenv('PDA_SITE'), Sys.getenv('PDA_SECRET'), 'digest'))
        obj<-jsonlite::fromJSON(file_path)
    }
}

    


#' @useDynLib PDA
#' @title PDA initialize
#' 
#' @usage pda_initialize <- function(mydata, broadcast=TRUE, control=pda_control)
#' @author Chongliang Luo, Steven Vitale
#' 
#' @param mydata local data in data frame
#' @param broadcast Logical, broadcast to the cloud? 
#' @param control PDA control
#'
#' @return  list(T_i = T_i, bhat_i = fit_i$coef, Vhat_i = summary(fit_i)$coef[,2]^2, site=control$mysite, site_size= nrow(mydata))
pda_initialize <- function(mydata,          
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
                 site = control$mysite,
                 site_size = nrow(mydata))
  }

   
  if(control$model=='ODACH'){
    # any diagnosis of heterogeneous baseline hazard?
    
    fit_i <- coxph(Surv(time, status) ~ ., data=mydata)
    
    init <- list(site = control$mysite,
                 site_size = nrow(mydata),
                 bhat_i = fit_i$coef,
                 Vhat_i = summary(fit_i)$coef[,2]^2   # cov matrix? vcov(fit_i)
                 )
  }
  
  
  if(control$model=='ODAL' | control$model=='ODALR'){ 
    fit_i <- glm(status ~ 0+., data=mydata,family = "binomial"(link = "logit"))
    init <- list(site = control$mysite,
                 site_size = nrow(mydata),
                 bhat_i = fit_i$coef,
                 Vhat_i = summary(fit_i)$coef[,2]^2)
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


#' @useDynLib PDA
#' @title PDA derivatives
#' 
#' @usage pda_derivatives(bbar, mydata, broadcast=TRUE, derivatives_ODAC_substep='first', control=pda_control)
#' @author Chongliang Luo, Steven Vitale
#' 
#' @param bbar  initial estimate
#' @param mydata local data in data frame
#' @param broadcast Logical, broadcast to the cloud? 
#' @param derivatives_ODAC_substep character, only for Cox regression, 'first' / 'second' indicate which substep of ODAC   
#' @param control PDA control
#' 
#' @details ONly for ODAC: step-2: calculate and broadcast 1st and 2nd order derivative at initial bbar
#'        for ODAC, this requires 2 substeps: 1st calculate summary stats (U, W, Z), 
#'        2nd calculate derivatives (logL_D1, logL_D2)
#'
#' @return  list(T_all=T_all, b_meta=b_meta, site=control$mysite, site_size = nrow(mydata), U=U, W=W, Z=Z, logL_D1=logL_D1, logL_D2=logL_D2)
pda_derivatives <- function(bbar = NULL,
                             mydata, 
                             broadcast = TRUE,
                             derivatives_ODAC_substep='first',   
                             control = pda_control){
  # data sanity check ...
  
  if(control$model=='ODAC'){
    px <- ncol(mydata) - 2
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
      derivatives <- list(T_all=T_all, b_meta=b_meta, site=control$mysite, site_size=nrow(mydata), U=U, W=W, Z=Z)
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
                          site=control$mysite, site_size = nrow(mydata),
                          logL_D1=logL_D1, logL_D2=logL_D2)
    }
  }
  
  
  if(control$model=='ODACH'){
    px <- ncol(mydata) - 2
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
    
    derivatives <- list(b_meta=sumstat_i$b_meta,  site=control$mysite, site_size = nrow(mydata),
                        logL_D1=logL_D1, logL_D2=logL_D2)

  }
  
  
  if(control$model=='ODAL' | control$model == "ODALR"){
    px <- ncol(mydata) - 1  # X includes intercept
    
    if(is.null(bbar)){ 
    # get b_meta as initial bbar
    bhat <- rep(0, px)
    vbhat <- rep(0, px)     # cov matrix?
    for(site_i in control$all_site){
      init_i <- pda_broadcast(obj_type= 'initialize',
                              upload=FALSE,
                              site_i=site_i, 
                              control=control) 
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
    }
    
    # 1st and 2nd derivatives
    status <- mydata$status
    X <- as.matrix(mydata[,-1])

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
    
    logL_D1 <- Lgradient(bbar,X,mydata$status)
    logL_D2 <- Lgradient2(bbar,X)
    
    derivatives <- list(
      site=control$mysite, 
      site_size = nrow(mydata),
      # b_init=bbar,
      logL_D1=logL_D1,
      logL_D2=logL_D2)
  }
  
  
  if(control$model=='ODALH'){
    
  }
  
  
  
  # broadcast to the cloud?
  if(broadcast){
    # if(not exist ...)
    if (control$model == "ODAL" | control$model == "ODALR"){obj_type <- 'derivatives'}
    # if(control$model=='ODAC' & derivatives_ODAC_substep == 'first')   
    else if(control$model=='ODAC' & derivatives_ODAC_substep == 'first'){obj_type <- 'derivatives_UWZ'}
    else{obj_type <- 'derivatives'}

    pda_broadcast(derivatives, obj_type, control=control)
  }else{
    cat('Derivatives not broadcasted. Please review the output and 
          use pda_broadcast() to manually upload to the cloud!')
  }
  
  return(derivatives)
}


#' @useDynLib PDA
#' @title PDA surrogate estimation
#' 
#' @usage pda_surrogate_est(bbar, mydata, broadcast=TRUE, control=pda_control)
#' @author Chongliang Luo, Steven Vitale
#' 
#' @param bbar  initial estimate
#' @param mydata local data in data frame
#' @param broadcast Logical, broadcast to the cloud? 
#' @param control PDA control
#' 
#' @details step-3: construct and solve surrogate logL at the master/lead site
#'
#' @return  list(btilde = sol$par, Htilde = sol$hessian, site=control$mysite, site_size=nrow(mydata))
pda_surrogate_est <- function(bbar = NULL,
                              mydata, 
                              broadcast = FALSE,        
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

    surr <- list(btilde = sol$par, Htilde = sol$hessian, site=control$mysite, site_size=nrow(mydata))
  }

  if(control$model=='ODAL' | control$model == "ODALR"){
    status <- mydata$status
    X <- as.matrix(mydata[,-1])
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
    if (control$model == "ODAL"){
      logL_all_D1 <- rep(0, px)
      logL_all_D2 <- matrix(0, px, px)
      N <- 0
      for(site_i in control$all_site){
        derivatives_i <- pda_broadcast(obj_type= 'derivatives',
                                       upload=FALSE,
                                       site_i=site_i,
                                       control=control)
        logL_all_D1 <- logL_all_D1 + derivatives_i$logL_D1*derivatives_i$site_size
        logL_all_D2 <- logL_all_D2 + derivatives_i$logL_D2*derivatives_i$site_size
        N <- N + derivatives_i$site_size
      }
      
      # initial beta
      if(is.null(bbar)) bbar <- control$beta_init  # derivatives_i$b_meta
      
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
      Y = mydata$status
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
      
    }else{
      logL_all_D1 <- rep(0, px)
      for(site_i in control$all_site){
        derivatives_i <- pda_broadcast(obj_type= 'derivatives',
                                       upload=FALSE,
                                       site_i=site_i,
                                       control=control)
        logL_all_D1 <- rbind(logL_all_D1,derivatives_i$logL_D1)
      }
      logL_all_D1 = logL_all_D1[-1,]
      
      # initial beta
      if(is.null(bbar)) bbar <- derivatives_i$b_meta
      
      #first order gradient
      Lgradient = function(beta,X,Y){
        design = X  # cbind(1,X)
        t(Y-expit(design%*%t(t(beta))))%*%design/length(Y)
      }
      #first-order surogate likelihood, suppose the local data are stored in Xlocal,Ylocal
      median_logL_all_D1 = apply(logL_all_D1, 2, median)
      logL_tilde = function(beta){
        - Lik(beta,X,mydata$status) - (median_logL_all_D1 - Lgradient(derivatives_i$b_meta,X,mydata$status))%*%beta
          # 2nd order?
      }
      
      # optimize the surrogate logL
      sol <- optim(par = bbar,
                   fn = logL_tilde,
                   # gr = logL_tilde_D1,
                   hessian = TRUE,
                   control = list(maxit=control$optim_maxit))
    }
    
    surr <- list(btilde = sol$par, Htilde = sol$hessian, site=control$mysite, site_size=nrow(mydata))
    ######################################################
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



#' @useDynLib PDA
#' @title PDA synthesize surrogate estimates from all sites, optional
#' 
#' @usage pda_synthesize(control=pda_control)
#' @author Chongliang Luo, Steven Vitale
#' 
#' @param control PDA control
#' 
#' @details Optional step-4: synthesize all the surrogate est btilde_i from each site, if step-3 from all sites is broadcasted
#'
#' @return  list(btilde=btilde,  Vtilde=Vtilde)
pda_synthesize <- function(control = pda_control){
  
  px <- length(control$risk_factor)
  K <- length(control$all_site)
  # if (control$model == "ODAL" | control$model == "ODALR"){
    # btilde_wt_sum <- rep(0, px+1)
    # wt_sum <- rep(0, px+1)     # cov matrix?
  # }else{
    btilde_wt_sum <- rep(0, px)
    wt_sum <- rep(0, px)     # cov matrix?
  # }
  
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
  
  cat("all surrogate estimates synthesized, no need to boradcast! \n")
  return(list(btilde=btilde, 
              Vtilde=Vtilde))
}



#' @useDynLib PDA
#' @import stats
#' @import Rcpp
#' @import survival 
#' @title PDA: Privacy-preserving Distributed Algorithm
#' 
#' @description  Fit Privacy-preserving Distributed Algorithms for linear, logistic, 
#'                Poisson and Cox PH regression with possible heterogeneous data across sites.
#' @usage pda(data = mydata, mysite = NULL)
#' @author Chongliang Luo, Steven Vitale, Jiayi Tong, Rui Duan, Yong Chen
#' 
#' @param mydata  Local IPD data in data frame, should include at least one column for the outcome and one column for the covariates 
#' @param mysite Character, site name as decided by all the sites when initialize the collaboration
#'
#'          
#' @references Jordan, Michael I., Jason D. Lee, and Yun Yang. "Communication-efficient distributed statistical inference." JASA (2018).
#' @examples
#'  ## Steve can you make an ODAl example here as we tested?
#' 
#' @export
pda <- function(data = mydata,
                mysite = Sys.getenv('PDA_SITE')){

  pda_control = pda_get('pda_control')
  cat('You are performing Privacy-preserving Distributed Algorithm (PDA, https://github.com/Penncil/PDA): \n')
  print(pda_control)
  # cat('your local analysis directory = ', mysite, '\n')
  if(is.null(mysite)) 
    stop(paste0('please specify your site name, one of ', paste(pda_control$all_site, collapse = ',')))
  else 
    cat('your site name = ', mysite, '\n')
  pda_control$mysite = mysite
  pda_control$model = pda_control$PDA_model
  
  n = nrow(data)
  formula<-as.formula(
    paste(pda_control$outcome,
    paste(pda_control$variables, collapse = " + "),
  sep = ' ~'))
  mf = model.frame(formula, data)
  if(pda_control$family=='cox'){  
    # myX = model.matrix(pda_control$formula, mf) 
    mydata = data.table(time=as.numeric(model.response(mf))[1:n], 
                        status=as.numeric(model.response(mf))[-c(1:n)], 
                        model.matrix(formula, mf)[,-1])
    pda_control$risk_factor = colnames(mydata)[-c(1:2)]
    # if(pda_control$heterogeneity==FALSE) pda_control$model = 'ODACH'
    # else pda_control$model = 'ODAC'
  }else{
    mydata = data.table(status=as.numeric(model.response(mf)), 
                        model.matrix(formula, mf))
    pda_control$risk_factor = colnames(mydata)[-1]
    # # family = 
    # # gaussian: DLM / DLMM, 
    # # binomial: ODAL / ODALH / ODALR, 
    # # poisson:  ODAP / ODAH / DPLR, 
    # # cox:      ODAC / ODACH
    # if(pda_control$heterogeneity==FALSE) pda_control$model = 'ODAL'
    # else pda_control$model = 'ODALH'
  }
  
  step = pda_control$step
  if(step==1 | step=='initialize'){
    output <- pda_initialize(mydata, control=pda_control)
    # print(output$bhat_i)
    # print(output$site)
    # print(output$site_size)
  } 
  
  # if(step==2 | step=='summary_stat'){
  #   output <- pda_summary_stat(bbar=NULL, mydata, control=control)
  #   print(output$b_meta) 
  # }
  
  if(step==2 | step=='derivatives') 
    output <- pda_derivatives(bbar=pda_control$beta_init, mydata, derivatives_ODAC_substep=derivatives_ODAC_substep, control=pda_control)
  
  if(step==3 | step=='surrogate_est'){
    output <- pda_surrogate_est(bbar=pda_control$beta_init, mydata, control=pda_control, broadcast = T)
    cat('\n', output$btilde)
    cat('\n Project accomplished, congratulations! \n')
    cat('If all sites uploaded their surrogate estimates, you can proceed to further synthesize them!')
    return(output)
  }
  
  if(step==4 | step=='synthesize'){
    output <- pda_synthesize(control=pda_control)
    print(output$btilde)
    return(output)
  } 
}


#' update the PDA control, used by the master site
#' @usage pda_control_update <- function(control_update=TRUE)
#' 
#' @param control_update Logical, update the PDA control?
#' 
#' @export
pda_control_update <- function(control_update=TRUE){
  pda_control = pda_get('pda_control')
  if(pda_control$step==1){
    init_i <- pda_broadcast(obj_type= 'initialize',
                            upload=FALSE,
                            site_i=pda_control$master_site, 
                            control=pda_control) 
    bhat <-init_i$bhat_i 
    vbhat <- init_i$Vhat_i
    for(site_i in pda_control$all_site){
      if(site_i!=pda_control$master_site){
        init_i <- pda_broadcast(obj_type= 'initialize',
                                upload=FALSE,
                                site_i=site_i, 
                                control=pda_control) 
        bhat = rbind(bhat, init_i$bhat_i)
        vbhat = rbind(vbhat, init_i$Vhat_i)
      }
    }
   
    #estimate from meta-analysis
    bmeta = apply(bhat/vbhat,2,function(x){sum(x, na.rm = T)})/apply(1/vbhat,2,function(x){sum(x, na.rm = T)})
    vmeta = 1/apply(1/vbhat,2,function(x){sum(x, na.rm = T)})
    res = list(bmeta = bmeta, vmeta = vmeta)
    cat('meta analysis (inverse variance weighted average) result:')
    print(res)
    pda_control$step = 2
    pda_control$beta_init = bmeta
    mes <- 'beta_init added, step=2 (derivatives)! \n'
  }else if(pda_control$step==2){
    res = ''
    pda_control$step = 3
    mes <- 'step=3 (surrogate_est)! \n'
  } else if(pda_control$step==3){
    res = ''
    pda_control$step = 4
    # synthesize surrogate est from each site?
    mes <- 'step=4 (synthesize)! \n'
  }
  
if(control_update==T){
  pda_broadcast(pda_control,
                obj_type= 'initialize',
                file_name = 'pda_control',
                upload=TRUE,
                # site_i=site_i, 
                control=pda_control)
  cat('pda_control has been updated on the cloud, ', mes)
  print(pda_control)
}

return(res)

}
 







################################# backup  ################################## 

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
