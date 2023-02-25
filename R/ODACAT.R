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

# set in pda()?
#NOTE: MLR models must contain intercept
ODACAT.steps <- c('initialize','derive','estimate','synthesize')
ODACAT.family <- 'multicategory'
# ODACAT.ordinal_categories<-FALSE

#' @useDynLib pda
#' @title ODACAT initialize
#' 
#' @usage ODACAT.initialize(ipdata,control,config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#' 
#' @return init
#' @keywords internal
ODACAT.initialize <- function(ipdata,control,config){
  n=nrow(ipdata)
  p=ncol(ipdata[,2:ncol(ipdata)])
  x= as.matrix(ipdata[,2:ncol(ipdata)])   # intercept will be added in model.fit... 
  y= as.matrix(ipdata[,1])
  
  #Proportional Odds Logistic Regression
  if(control$ordinal_categories==FALSE){
    b.t=rep(0,(p+1)*(control$number_outcome_categories-1))
    fit.p=model.fit(x=x, y=y, b0=b.t, model = 'mlr')
    init=list(site = config$site_id,
              site_size = n,
              bhat_i = as.vector(fit.p$res$par),
              Vhat_i = diag(solve(fit.p$res$hessian*n)),
              beta=NA,  
              zeta=NA, 
              theta=NA) 
  }else{
    #Parameters are estimated in terms of beta and theta
    b.t=rep(0,p+(control$number_outcome_categories-1))
    fit.p=model.fit(x=x, y=y, b0=b.t, model = 'polr')
    init=list(site = config$site_id,
              site_size = n,
              bhat_i = as.vector(fit.p$res$par),
              Vhat_i = diag(solve(fit.p$res$hessian*n)),
              beta=fit.p$beta,
              zeta=fit.p$zeta,   # transform back to int terms, has order restriction
              theta=fit.p$theta) # no order restriction
  }
  return(init)
}

#' @useDynLib pda
#' @title ODACAT derivatives
#' 
#' @usage ODACAT.derive(ipdata,control,config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#'
#' @return  list(site=config$site_id, site_size = nrow(ipdata), logL_D1=logL_D1, logL_D2=logL_D2)
#' @keywords internal
ODACAT.derive <- function(ipdata,control,config){
  n=nrow(ipdata)  
  p=ncol(ipdata[,2:ncol(ipdata)])
  if(control$ordinal_categories==FALSE){
    pp=(p+1)*(control$number_outcome_categories-1)
  }else{
    pp=p+(control$number_outcome_categories-1)
  }
  
  bhat <- control$bhat
  vbhat <- control$vbhat
  #Meta Analysis
  bhat <- rep(0, pp)
  vbhat <- rep(0, pp)     
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
  
  # meta_res=meta.fit(ipdata,control,config)  
  # bbar=meta_res$bbar #b_meta
  # vbar=meta_res$vmeta
  
  x= as.matrix(ipdata[,2:ncol(ipdata)])   # intercept will be added in model.fit... 
  y= as.matrix(ipdata[,1])
  
  #First and Second Order Gradients
  if(control$ordinal_categories==FALSE){
    gradients=model.logL.d(g=bbar, x=x, y=y, model="mlr", grad=T, hess=T)
  }else{
    gradients=model.logL.d(g=bbar, x=x, y=y, model="polr", grad=T, hess=T)
  }
  
  derivatives <- list(
    site=config$site_id, 
    site_size = n,
    logL_D1=gradients$logL.grad,
    logL_D2=gradients$logL.hess)
  
  return(derivatives)
}


#' @useDynLib pda
#' @title PDA surrogate estimation
#' 
#' @usage ODACAT.estimate(ipdata,control,config)
#' @param ipdata local data in data frame
#' @param control PDA control
#' @param config cloud configuration
#' 
#' @details step-3: construct and solve surrogate logL at the master/lead site
#' @return  list(btilde = sol$par, Htilde = sol$hessian, site=control$mysite, site_size=nrow(ipdata))
#' @keywords internal
ODACAT.estimate <- function(ipdata,control,config) {
  
  n_site=nrow(ipdata)
  p=ncol(ipdata[,2:ncol(ipdata)])
  K=length(control$sites)
  
  if(control$ordinal_categories==FALSE){
    pn=(p+1)*(control$number_outcome_categories-1)
  }else{
    pn=p+(control$number_outcome_categories-1)
  }
  
  logL_N_D1_beta_bar <- matrix(NA, K, pn)
  logL_N_D2_beta_bar <- array(NA, c(K, pn, pn))
  nn=rep(0,K)
  for(i in 1:length(control$sites)){   
    site_i=control$sites[i]
    derivatives_i=pdaGet(paste0(site_i,'_derive'),config)
    logL_N_D1_beta_bar[i,]=derivatives_i$logL_D1
    logL_N_D2_beta_bar[i,,]=derivatives_i$logL_D2
    nn[i]=derivatives_i$site_size
  }
  n=sum(nn)
  
  # initial beta
  bbar=control$beta_init  # derivatives_i$b_meta
  
  #Local Site
  local_site=config$site_id
  local_site_index=match(local_site,control$sites)
  
  
  #Surrogate Likelihood Estimation
  x= as.matrix(ipdata[,2:ncol(ipdata)])   # intercept will be added in model.fit... 
  y= as.matrix(ipdata[,1])
  if(control$ordinal_categories==FALSE){
    surr_res=odacat.fit(x=x, y=y, #mydata=ipdata,
                        id.local=local_site_index,beta_bar=bbar,logL_N_D1_beta_bar=logL_N_D1_beta_bar,
                        logL_N_D2_beta_bar=logL_N_D2_beta_bar,nn=nn, model='mlr')
    vcov=solve(surr_res$sol.odacat2.hess*n)
    surr <- list(btilde = surr_res$sol.odacat2$par, 
                 btilde.se = sqrt(diag(vcov)), 
                 btilde.theta=NULL,
                 btilde.theta.se=NULL,
                 site=config$site_id, 
                 site_size=n_site,
                 vcov=vcov)
  }else{
    surr_res=odacat.fit(x=x, y=y, # mydata=ipdata,
                        id.local=local_site_index,beta_bar=bbar,logL_N_D1_beta_bar=logL_N_D1_beta_bar,
                        logL_N_D2_beta_bar=logL_N_D2_beta_bar,nn=nn, model='polr')
    vcov.betatheta=solve(surr_res$sol.odacat2.hess*n)
    surr_res_param_repar=repar2(g=surr_res$sol.odacat2$par, g.covar=vcov.betatheta, px=p)
    surr <- list(btilde = c(surr_res_param_repar$beta,surr_res_param_repar$zeta), 
                 btilde.se = sqrt(c(surr_res_param_repar$beta.var,surr_res_param_repar$zeta.var)),
                 btilde.theta=surr_res_param_repar$theta,
                 btilde.theta.se=sqrt(surr_res_param_repar$theta.var),
                 site=config$site_id, 
                 site_size=n_site,
                 vcov=surr_res_param_repar$vcov)
  }
  
  return(surr)
}

#' @useDynLib pda
#' @title PDA synthesize surrogate estimates from all sites, optional
#' 
#' @usage ODACAT.synthesize(ipdata,control,config)
#' @param ipdata local data in data frame
#' @param control pda control
#' @param config pda cloud configuration
#' 
#' @details Optional step-4: synthesize all the surrogate est btilde from each site, if step-3 from all sites is broadcasted
#'
#' @return  list(btilde=btilde, Vtilde=Vtilde)
#' @keywords internal
ODACAT.synthesize <- function(ipdata,control,config) {
  p=ncol(ipdata[,2:ncol(ipdata)])
  if(control$ordinal_categories==FALSE){
    px=(p+1)*(control$number_outcome_categories-1)
  }else{
    px=p+(control$number_outcome_categories-1)
  }
  
  K <- length(control$sites)
  btilde_matrix=matrix(0,K,px)
  vcov_sum <- rep(0, px)  
  btilde.var_matrix=matrix(0,K,px)
  
  
  
  for(site_i in control$sites){
    surr_i <- pdaGet(paste0(site_i,'_estimate'),config)
    btilde_matrix[i,]=surr_i$btilde
    btilde.var_matrix[i,]=diag(surr_i$vcov)
    vcov_sum <- vcov_sum + surr_i$vcov
  }
  
  # inv-Var weighted average est, and final Var = average Var-tilde
  btilde=apply((1/btilde.var_matrix) *  btilde_matrix, 2,sum)/apply((1/btilde.var_matrix), 2, sum)
  Vtilde=vcov_sum/K
  
  
  message("all surrogate estimates synthesized, no need to broadcast! ")
  return(list(btilde = btilde,
              Vtilde = Vtilde))
}


