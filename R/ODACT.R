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

# Rcpp::sourceCpp('pda/src/odact.cpp')
ODACT.steps<-c('initialize','derive','estimate','synthesize')  # ODACH # 'deriveUWZ',
ODACT.family<-'cox'



#' @useDynLib pda
#' @title ODACT initialize
#' 
#' @usage ODACT.initialize(ipdata, control, config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#' 
#' @references  Liang CJ, Luo C, Kranzler HR, Bian J, Chen Y. Communication-efficient federated learning of temporal effects on opioid use disorder with data from distributed research networks. J Am Med Inform Assoc. 2025 Apr 1;32(4):656-664. doi: 10.1093/jamia/ocae313. PMID: 39864407; PMCID: PMC12005629. 
#' @return  list(bhat_i, Vhat_i, site, site_size)
#' @keywords internal
ODACT.initialize <- function(ipdata,control,config){
  # if(control$heterogeneity == FALSE){
  #   T_i <- sort(unique(ipdata$time[ipdata$status==TRUE]))
  # }else{
  #   T_i <- NA
  # }
  
  # handle data degeneration (e.g. missing categories in some site). This could be in pda()?
  px = ncol(ipdata) - 2 
  h = control$bandwidth
  evalt = control$times
  nT = length(evalt)
  col_deg = apply(ipdata[,-c(1:2)],2,var)==0    # degenerated X columns...
  ipdata_i = ipdata[,-(which(col_deg)+2),with=F] 
  
  # Cox with beta(t): local constant partial likelihood estimate
  fit_i <- seq_fit_list(data.frame(ipdata_i), fn=llpl, times=evalt, h=h, betabar=rep(0,px), hessian=T) 
  
  if(!is.null(fit_i)){
    # for degenerated X, coef=0, var=Inf
    bhat_i = matrix(0,px,nT)    # px * nT 
    Vhat_i = matrix(Inf,px,nT)
    bhat_i[!col_deg,] <- sapply(fit_i, function(a) a$par) 
    Vhat_i[!col_deg,] <- sapply(fit_i, function(a) diag(solve(a$hessian+diag(1e-7,px))*0.6/h) )  
    
    init <- list(#T_i = T_i,
                 bhat_i = bhat_i,
                 Vhat_i = Vhat_i,   #   
                 site = config$site_id,
                 site_size = nrow(ipdata))
    init$Vhat_i[init$Vhat_i==0] = NA
  } else{
    init <- list(#T_i = T_i,
                 bhat_i = NA,
                 Vhat_i = NA,   
                 site = config$site_id,
                 site_size = nrow(ipdata))
  }
  
  return(init)
}


#' @useDynLib pda
#' @title Generate pda ODACT derivatives
#' 
#' @usage ODACT.derive(ipdata, control, config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#' 
#' @details Calculate and broadcast 1st and 2nd order derivative at initial bbar for ODACT 
#'
#' @import Rcpp  
#' @return  list(b_meta=b_meta, site=control$mysite, site_size = nrow(ipdata), logL_D1=logL_D1, logL_D2=logL_D2)
#' @keywords internal
ODACT.derive <- function(ipdata,control,config){
  px <- ncol(ipdata) - 2  
  time <- ipdata$time
  status <- ipdata$status
  X <- as.matrix(ipdata[,-c(1,2)])
  n <- length(time)
  bbar = control$beta_init # px * nT
  h = control$bandwidth
  evalt = control$times

  # lists (nT) of gradient/hessian...
  logL_D1 <- lapply(as.list(evalt), function(x) llplg(bbar[,findInterval(x, evalt)], time, status, X, tt=x, h=h, order.time=F))     
  logL_D2 <- lapply(as.list(evalt), function(x) llplh(bbar[,findInterval(x, evalt)], time, status, X, tt=x, h=h, order.time=F))
    
  derivatives <- list(b_init=bbar, 
                      site=config$site_id, site_size = nrow(ipdata),
                      logL_D1=logL_D1, logL_D2=logL_D2) 
  
  return(derivatives)
}


#' @useDynLib pda
#' @title PDA ODACT surrogate estimation
#' 
#' @usage ODACT.estimate(ipdata, control, config)
#' @param ipdata local data in data frame
#' @param control pda control
#' @param config cloud config
#' @import data.table
#' 
#' @details step-3: construct and solve surrogate logL at the Lead site
#' @import Rcpp  
#' @return  list(btilde = sol$par, Htilde = sol$hessian, site=control$mysite, site_size=nrow(ipdata))
#' @keywords internal
ODACT.estimate <- function(ipdata,control,config) {
  # data sanity check ...
  time <- ipdata$time
  status <- ipdata$status
  X <- as.matrix(ipdata[,-c(1,2)])
  n <- length(time)
  px <- ncol(X) 
  h = control$bandwidth
  evalt = control$times
  nT = length(evalt)
  K = length(control$sites)
  
  # download derivatives of other sites from the cloud
  # calculate 2nd order approx of the total logL  
  logL_all_D1 <- matrix(0, nT, px)
  logL_all_D2 <- array(0, c(nT, px, px))
  N <- 0
  for(i in 1:K){
    derivatives_i <- pdaGet(paste0(control$sites[i],'_derive'),config)
    logL_all_D1 <- logL_all_D1 + derivatives_i$logL_D1 
    logL_all_D2 <- logL_all_D2 + derivatives_i$logL_D2 
    N <- N + derivatives_i$site_size
  } 
  
  # initial beta
  bbar <- control$beta_init
  btilde = setilde = matrix(0, px, nT)
  for(it in 1:nT){ # at each time point
    tt = evalt[it]
    # logL at local site 
    logL_local <- function(b) llpl(b, time, status, X, tt, h=h, order.time=F)
    logL_local_D1 <- function(b) llplg(b, time, status, X, tt, h=h, order.time=F)     
    logL_local_D2 <- function(b) llplh(b, time, status, X, tt, h=h, order.time=F)
    
    # surrogate log-L and its gradient
    logL_diff_D1 <- logL_all_D1[it,] / N - logL_local_D1(bbar[,it]) / n
    logL_diff_D2 <- logL_all_D2[it,,] / N - logL_local_D2(bbar[,it]) / n
    logL_tilde <- function(b) c(logL_local(b) / n + sum(b * logL_diff_D1) + 1/2 * t(b-bbar[,it]) %*% logL_diff_D2 %*% (b-bbar[,it]))

    # optimize the surrogate logL 
    sol <- optim(par = bbar[,it], 
                 fn = logL_tilde,
                 # gr = logL_tilde_D1,
                 hessian = TRUE,
                 method = control$optim_method,
                 control = list(maxit=control$optim_maxit))
    btilde[,it] = sol$par
    # var estimate: by inv hessian 
    setilde[,it] = sqrt(diag(solve(sol$hessian+diag(1e-7,px)))/N)
  }
  
  surr <- list(btilde = btilde, setilde=setilde,  
               site=config$site_id, site_size=nrow(ipdata) )
  return(surr)
}



#' @useDynLib pda
#' @title PDA synthesize surrogate estimates from all sites, optional
#' 
#' @usage ODACT.synthesize(ipdata, control, config)
#' @param ipdata local data in data frame
#' @param control pda control
#' @param config cloud config
#' 
#' @details Optional step-4: synthesize all the surrogate est btilde_i from each site, if step-3 from all sites is broadcasted
#' @import Rcpp  
#' @return  list(btilde=btilde,  Vtilde=Vtilde)
#' @keywords internal
ODACT.synthesize <- function(ipdata,control,config) { 
  n <- length(time) 
  h = control$bandwidth
  evalt = control$times
  nT = length(evalt)
  px <- length(control$risk_factor)
  K <- length(control$sites)
  btilde = array(0, c(px, nT, K))
  setilde = array(0, c(px, nT, K))
  
  for(k in 1:K){
    surr_i <- pdaGet(paste0(control$sites[k],'_estimate'),config)
    btilde[,,k] = surr_i$btilde 
    setilde[,,k] = surr_i$setilde 
  }
  b_wt_sum <- apply(btilde / (setilde^2), c(1,2), sum, na.rm=T)  
  wt_sum <- apply(1 / (setilde^2), c(1,2), sum, na.rm=T) 
  
  # inv-Var weighted average est, and final Var = average Var-tilde
  btilde <- b_wt_sum / wt_sum  
  setilde <- sqrt(K / wt_sum)  
  
  message("all surrogate estimates synthesized, no need to broadcast! ")
  return(list(btilde=btilde, setilde=setilde ))
}
