# https://style.tidyverse.org/functions.html#naming

# require(survival)
# require(data.table)
# Rcpp::sourceCpp('src/rcpp_coxph.cpp')
## broadcast: upload/download shared info to/from the cloud folder

#' @useDynLib pda
#' @title Function to upload object to cloud
#' 
#' @usage pda_put(obj,name)
#' @author Chongliang Luo, Steven Vitale
#' 
#' @param obj R object to upload
#' @param name of file
#'
#' @return  
pda_put <- function(obj,name,config){
    obj_Json <- jsonlite::toJSON(obj)
    file_name <- paste0(name, '.json')
    print(paste("Put",file_name,"on public cloud:"))
    print(obj_Json)
#    if(interactive()) {
#      authorize = menu(c("Yes", "No"), title="Allow upload?")
#    } else {
      authorize = "1"
#    }
    if (authorize != 1) {
      print("file not uploaded.")
      return(FALSE)
    }
    # the file to upload
    if (is.character(config$dir)) {
        file_path <- paste0(config$dir,'/', file_name)
    } else {
        file_path <- paste0(tempdir(),'/', file_name)
    }
    write(obj_Json, file_path)
    print(paste("wrote to",file_path))
    if (is.character(config$pda_uri)) {
        # create the url target of the file
        url <- file.path(config$uri, file_name)
        # webdav PUT request to send a file to cloud
        r<-httr::PUT(url, body = httr::upload_file(file_path), httr::authenticate(config$site_id, config$site_secret, 'digest'))
        print("file uploaded")
    }
}


#' @useDynLib pda
#' @title Function to download json and return as object)
#' 
#' @usage pda_get(name)
#' @author Chongliang Luo, Steven Vitale
#' 
#' @param name of file
#' @return  
#' @export
pda_get <- function(name,config){
    file_name <- paste0(name, '.json')
    print(paste("Get",file_name,"on public cloud:"))
    # the file to upload
    if (is.character(config$dir)) {
        file_path <- paste0(config$dir,'/', file_name)
    } else {
        file_path <- paste0(tempdir(),'/', file_name)
    }
    if (is.character(config$uri)) {
        url <- file.path(config$uri, file_name)
        #write the file from GET request to file_path
        res<-httr::GET(url, httr::write_disk(file_path, overwrite = TRUE), httr::authenticate(config$site_id, config$site_secret, 'digest'))
    } 
    obj<-jsonlite::fromJSON(file_path)
    return(obj)
}

    


#' @useDynLib pda
#' @title PDA initialize
#' 
#' @usage pda_initialize <- function(ipdata, broadcast=TRUE, control=control)
#' @author Chongliang Luo, Steven Vitale
#' @param ipdata local data in data frame
#' @param control pda control
#' @return  list(T_i = T_i, bhat_i = fit_i$coef, Vhat_i = summary(fit_i)$coef[,2]^2, site=control$mysite, site_size= nrow(ipdata))
pda_initialize <- function(ipdata,control,config){
  # data sanity check ...
  
  # if(!any(names(ipdata)[1:2] == c('time', 'status')))
  #   error('ipdata columns should be (time, status, covariates)')
  # if(!any(is.numeric(ipdata)))
  #   error('ipdata need to be numeric, please create dummy variables if necessary')
  
  if(control$model=='ODAC'){
    T_i <- sort(unique(ipdata$time[ipdata$status==TRUE]))
    fit_i <- coxph(Surv(time, status) ~ ., data=ipdata)
    
    init <- list(T_i = T_i,
                 bhat_i = fit_i$coef,
                 Vhat_i = summary(fit_i)$coef[,2]^2,   # cov matrix? vcov(fit_i)
                 site = config$site_id,
                 site_size = nrow(ipdata))
  }

   
  if(control$model=='ODACH'){
    # any diagnosis of heterogeneous baseline hazard?
    
    fit_i <- coxph(Surv(time, status) ~ ., data=ipdata)
    
    init <- list(site = config$site_id,
                 site_size = nrow(ipdata),
                 bhat_i = fit_i$coef,
                 Vhat_i = summary(fit_i)$coef[,2]^2   # cov matrix? vcov(fit_i)
                 )
  }
  
  
  if(control$model=='ODAL' | control$model=='ODALR'){ 
    fit_i <- glm(status ~ 0+., data=ipdata,family = "binomial"(link = "logit"))
    init <- list(site = config$site_id,
                 site_size = nrow(ipdata),
                 bhat_i = fit_i$coef,
                 Vhat_i = summary(fit_i)$coef[,2]^2)
  }
  
  
  if(control$model=='ODALH'){
    
  }
  
  
  #  pda_put(init, paste_0(config$site_id,'_initialize',config))
  return(init)
}


#' @useDynLib pda
#' @title PDA derivatives
#' 
#' @usage pda_derivatives(bbar, ipdata, broadcast=TRUE, derivatives_ODAC_substep='first', control=control)
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
pda_derivatives <- function(ipdata,control,config){
  # data sanity check ...
  if(control$model=='ODAL' | control$model == "ODALR"){
    px <- ncol(ipdata) - 1  # X includes intercept
    # get b_meta as initial bbar
    bhat <- rep(0, px)
    vbhat <- rep(0, px)     # cov matrix?
    for(site_i in control$all_site){
      init_i <- pda_get(paste0(site_i,'_initialize'),config)
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
pda_surrogate_est <- function(ipdata,control,config) {
  # data sanity check ...
  if(control$model=='ODAC' | control$model=='ODACH'){
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
    for(site_i in control$all_site){
      derivatives_i <- pda_get(paste0(site_i,'_derivatives'),config)
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

    surr <- list(btilde = sol$par, Htilde = sol$hessian, site=config$site_id, site_size=nrow(ipdata))
  }

  if(control$model=='ODAL' | control$model == "ODALR"){
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
    if (control$model == "ODAL"){
      logL_all_D1 <- rep(0, px)
      logL_all_D2 <- matrix(0, px, px)
      N <- 0
      for(site_i in control$all_site){
        derivatives_i <- pda_get(paste0(site_i,'_derivatives'),config)
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
      
    }else{
      logL_all_D1 <- rep(0, px)
      for(site_i in control$all_site){
        derivatives_i <- pda_get(paste0(site_i,'_derivatives'),config)
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
        - Lik(beta,X,ipdata$status) - (median_logL_all_D1 - Lgradient(derivatives_i$b_meta,X,ipdata$status))%*%beta
          # 2nd order?
      }
      
      # optimize the surrogate logL
      sol <- optim(par = bbar,
                   fn = logL_tilde,
                   # gr = logL_tilde_D1,
                   hessian = TRUE,
                   control = list(maxit=control$optim_maxit))
    }
    
    surr <- list(btilde = sol$par, Htilde = sol$hessian, site=config$site_id, site_size=nrow(ipdata))
    ######################################################
  }

    
  if(control$model=='ODALH'){
  
  }
  
  # broadcast to the cloud?
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
pda_synthesize <- function(ipdata,control,config) {
  
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
    surr_i <- pda_get(paste0(site_i,'_derivatives'),config)
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



#' @useDynLib pda
#' @import stats
#' @import Rcpp
#' @import survival 
#' @title PDA: Privacy-preserving Distributed Algorithm
#' 
#' @description  Fit Privacy-preserving Distributed Algorithms for linear, logistic, 
#'                Poisson and Cox PH regression with possible heterogeneous data across sites.
#' @usage pda(data = ipdata, mysite = NULL)
#' @author Chongliang Luo, Steven Vitale, Jiayi Tong, Rui Duan, Yong Chen
#' 
#' @param ipdata  Local IPD data in data frame, should include at least one column for the outcome and one column for the covariates 
#' @param mysite Character, site name as decided by all the sites when initialize the collaboration
#'
#'          
#' @references Jordan, Michael I., Jason D. Lee, and Yun Yang. "Communication-efficient distributed statistical inference." JASA (2018).
#' @examples
#'  ## Steve can you make an ODAl example here as we tested?
#' 
#' @export
pda <- function(ipdata,config){
  control = pda_get('control',config)
  cat('You are performing Privacy-preserving Distributed Algorithm (PDA, https://github.com/Penncil/pda): \n')
  cat('your site = ', config$site_id, '\n')
  n = nrow(ipdata)
  formula<-as.formula(
    paste(control$outcome,
    paste(control$variables, collapse = " + "),
  sep = ' ~'))
  mf = model.frame(formula, ipdata)
  if(control$family=='cox'){  
    ipdata = data.table(time=as.numeric(model.response(mf))[1:n], 
                        status=as.numeric(model.response(mf))[-c(1:n)], 
                        model.matrix(formula, mf)[,-1])
    control$risk_factor = colnames(ipdata)[-c(1:2)]
  }else{
    ipdata = data.table(status=as.numeric(model.response(mf)), 
                        model.matrix(formula, mf))
    control$risk_factor = colnames(ipdata)[-1]
    # # family = 
    # # gaussian: DLM / DLMM, 
    # # binomial: ODAL / ODALH / ODALR, 
    # # poisson:  ODAP / ODAH / DPLR, 
    # # cox:      ODAC / ODACH
    # if(control$heterogeneity==FALSE) control$model = 'ODAL'
    # else control$model = 'ODALH'
  }
  if(control$step=='initialize'){
    return(pda_initialize(ipdata, control, config))
  } else if(control$step=='derivatives') {
    return(pda_derivatives(ipdata, control, config))
  } else if(control$step=='derivatives_UWZ') {
    return(pda_derivatives_UWZ(ipdata, control, config))
  } else if(control$step=='derivatives') {
    return(pda_derivatives(ipdata, control, config))
  } else if(control$step=='surrogate_est'){
    return(pda_surrogate_est(ipdata,control,config))
  } else if(control$step=='synthesize'){
    return(pda_synthesize(ipdata,control,config))
  } 
}


#' update the PDA control, used by the master site
#' @usage control_update <- function(control_update=TRUE)
#' 
#' @param control_update Logical, update the PDA control?
#' 
#' @export
control_update <- function(config){
  control = pda_get('control',config)
  if(control$step=="initialize"){
    init_i <- pda_get(paste0(control$lead_site,'_initialize'),config)
    bhat <-init_i$bhat_i 
    vbhat <- init_i$Vhat_i
    for(site_i in control$all_site){
      if(site_i!=control$lead_site){
        init_i <- pda_get(paste0(site_i,'_initialize'),config)
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
    control$step = "derivatives"
    control$beta_init = bmeta
    mes <- 'beta_init added, step=2 (derivatives)! \n'
  }else if(control$step=="derivatives"){
    res = ''
    control$step = "surrogate_est"
    mes <- 'step=3 (surrogate_est)! \n'
  } else if(control$step=="surrogate_est"){
    res = ''
    control$step = "synthesize" 
    # synthesize surrogate est from each site?
    mes <- 'step=4 (synthesize)! \n'
  }
  
  cat(mes)
  return(control)

}
 







################################# backup  ################################## 

pda_main <- function(ipdata = ipdata,
                     step = 1,                       # c('initialize', 'derivatives', 'surrogate_est', 'synthesize'),
                     derivatives_ODAC_substep=NULL,  # c('first', 'second'),  only for control$model=='ODAC' and step==2
                     control = control){
  if(step=='initialize'){
    output <- pda_initialize(ipdata, control=control)
    print(output$bhat_i)
    print(output$site)
    print(output$site_size)
  } 
  
  # if(step==2 | step=='summary_stat'){
  #   output <- pda_summary_stat(bbar=NULL, ipdata, control=control)
  #   print(output$b_meta) 
  # }
  
  if(step=='derivatives') 
    output <- pda_derivatives(bbar=NULL, ipdata, derivatives_ODAC_substep=derivatives_ODAC_substep, control=control)
  
  
  if(step=='surrogate_est'){
    output <- pda_surrogate_est(bbar=NULL, ipdata, control=control)
    cat('\n', output$btilde)
    return(output)
  }
  
  if(step=='synthesize'){
    output <- pda_synthesize(control=control)
    print(output$btilde)
    return(output)
  }
  
  # return(output)
}
