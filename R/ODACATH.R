# Copyright 2023 Penn Computing Inference Learning (PennCIL) lab
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
ODACATH.steps <- c('initialize','derive','estimate','synthesize')
ODACATH.family <- 'multicategory'
# ODACAT.ordinal_categories<-FALSE

#' @useDynLib pda
#' @title ODACATH initialize
#' 
#' @usage ODACATH.initialize(ipdata,control,config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#' 
#' @return init
#' @keywords internal
ODACATH.initialize <- function(ipdata,control,config){
  #Estimate coefficients for predictors of interest and site specific intercepts
  n=nrow(ipdata)
  p=ncol(ipdata[,2:ncol(ipdata)])
  x= as.matrix(ipdata[,2:ncol(ipdata)])   
  y= as.matrix(ipdata[,1])
  
  
  #Proportional Odds Logistic Regression
  if(control$ordinal_categories==FALSE){ #MLR
    l_tot=(p+1)*(control$number_outcome_categories-1)
    beta_length=p*(control$number_outcome_categories-1)
    eta_length=control$number_outcome_categories-1
    beta0_record = rep(0,beta_length) #Shared beta coefficients by site
    eta0_record = rep(0,eta_length) #Site specific eta coefficients by site
    eta_indices=sapply(1:(control$number_outcome_categories-1),function(l){1+(l-1)*(p+1)})
    beta_indices=seq(1:l_tot)
    beta_indices=beta_indices[-eta_indices]
    
    sl1 = function(g) {Lik2(g=g,x=x,y=y,model="mlr")} #Log Likelihood
    theta1_0=optim(rep(0,times=l_tot), sl1,method="BFGS", hessian=T, control = list(maxit=500,reltol=1e-8)) #Estimation of parameters 
    beta0_record=theta1_0$par[beta_indices] #Recording estimated shared beta coefficients
    eta0_record=theta1_0$par[eta_indices] #Recording estimated site specific beta coefficients
    
    ### compute the weight (Can alternatively use the inverse hessian of the results from nlminb)
    ###Variance of coefficients computed as \sigma^2=(X'WX) where W=diag(p(1-p)) and p=exp(xb)/(1+exp(xb))
    Vall=diag(solve(theta1_0$hessian*n))
    V_beta=Vall[beta_indices]
    V_eta=Vall[eta_indices]
    
    init=list(site = config$site_id,
              site_size = n,
              bhat_i =beta0_record ,
              Vhat_i =V_beta ,
              bhat_eta_i=eta0_record,
              Vhat_eta_i=V_eta,
              beta=NA,  
              zeta=NA, 
              theta=NA) 
  }else{ #POLR
    l_tot=p+(control$number_outcome_categories-1)
    beta_length=p
    eta_length=control$number_outcome_categories-1
    beta0_record = rep(0,beta_length) #Shared beta coefficients by site
    eta0_record = rep(0,eta_length) #Site specific eta coefficients by site
    eta_indices=(1:(control$number_outcome_categories-1))+p
    beta_indices=seq(1:l_tot)
    beta_indices=beta_indices[-eta_indices]
    
    sl1 = function(g) {Lik2(g=g,x=x,y=y,model="polr")} #Log Likelihood
    theta1_0=optim(rep(0,times=l_tot), sl1,method="BFGS", hessian=T, control = list(maxit=500,reltol=1e-8)) #Estimation of parameters 
    beta0_record=theta1_0$par[beta_indices] #Recording estimated shared beta coefficients
    eta0_record=theta1_0$par[eta_indices] #Recording estimated site specific beta coefficients
    
    ### compute the weight (Can alternatively use the inverse hessian of the results from nlminb)
    ###Variance of coefficients computed as \sigma^2=(X'WX) where W=diag(p(1-p)) and p=exp(xb)/(1+exp(xb))
    Vall=diag(solve(theta1_0$hessian*n))
    V_beta=Vall[beta_indices]
    V_eta=Vall[eta_indices]
    
    tt=repar(g=theta1_0$par, g.covar=solve(theta1_0$hessian*n), px=p)
    
    init=list(site = config$site_id,
              site_size = n,
              bhat_beta_i =beta0_record ,
              Vhat_beta_i =V_beta ,
              bhat_eta_i=eta0_record,
              Vhat_eta_i=V_eta,
              beta=tt$beta,  
              zeta=tt$zeta, 
              theta=tt$theta)
  }
  return(init)
}

#' @useDynLib pda
#' @title ODACATH derivatives
#' 
#' @usage ODACATH.derive(ipdata,control,config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#'
#' @return  list(site=config$site_id, site_size = n, S_site=S_site, eta=eta_mat[site,])
#' @keywords internal
ODACATH.derive <- function(ipdata,control,config){
  n=nrow(ipdata)
  p=ncol(ipdata[,2:ncol(ipdata)])
  x= as.matrix(ipdata[,2:ncol(ipdata)])   
  y= as.matrix(ipdata[,1])
  if(control$ordinal_categories==FALSE){
    model = "mlr"
    l_tot=(p+1)*(control$number_outcome_categories-1)
    beta_length=p*(control$number_outcome_categories-1)
    eta_length=control$number_outcome_categories-1
    beta0_record = rep(0,beta_length) #Shared beta coefficients by site
    eta0_record = rep(0,eta_length) #Site specific eta coefficients by site
    eta_indices=sapply(1:(control$number_outcome_categories-1),function(l){1+(l-1)*(p+1)})
    # beta_indices=seq(1:l_tot)
    # beta_indices=beta_indices[-eta_indices]
  }else{
    model = "polr"
    l_tot=p+(control$number_outcome_categories-1)
    beta_length=p
    eta_length=control$number_outcome_categories-1
    beta0_record = rep(0,beta_length) #Shared beta coefficients by site
    eta0_record = rep(0,eta_length) #Site specific eta coefficients by site
    eta_indices=(1:(control$number_outcome_categories-1))+p
    # beta_indices=seq(1:l_tot)
    # beta_indices=beta_indices[-eta_indices]
  }
  
 
  bbar = control$beta_init
  eta_mat = matrix(unlist(control$bhat_eta),ncol = length(eta_indices), byrow = TRUE)
  
  x= as.matrix(ipdata[,2:ncol(ipdata)])   # intercept will be added in model.fit... 
  y= as.matrix(ipdata[,1])
  
  #Efficient Score within Site evaluated at beta-bar
  #Local Site
  site_name=config$site_id
  site=match(site_name,control$sites) #Site index
  model = control$model
  #site =as.numeric(substr(siteid,nchar(siteid), nchar(siteid)))
  S_site=regular_es(beta=bbar,eta=eta_mat[site,],x=x,y=y,eta_indices=eta_indices,model=model)*n #Will divide the sum across all sites for N
  
  derivatives <- list(
    site=config$site_id, 
    site_size = n,
    S_site=S_site,
    eta=eta_mat[site,],
    bbar=bbar)
  
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
#' @details step-3: construct and solve surrogate efficient score at the master/lead site
#' @return  list(btilde=betanew, btilde.se=beta_SE,eta_mat=eta_mat,eta_mat_theta=NULL,site=config$site_id, site_size=n_site)
#' @keywords internal
ODACATH.estimate <- function(ipdata,control,config) {
  
  n=nrow(ipdata)
  p=ncol(ipdata[,2:ncol(ipdata)])
  x= as.matrix(ipdata[,2:ncol(ipdata)])   
  y= as.matrix(ipdata[,1])
  if(control$ordinal_categories==FALSE){
    l_tot=(p+1)*(control$number_outcome_categories-1)
    beta_length=p*(control$number_outcome_categories-1)
    eta_length=control$number_outcome_categories-1
    beta0_record = rep(0,beta_length) #Shared beta coefficients by site
    eta0_record = rep(0,eta_length) #Site specific eta coefficients by site
    eta_indices=sapply(1:(control$number_outcome_categories-1),function(l){1+(l-1)*(p+1)})
    beta_indices=seq(1:l_tot)
    beta_indices=beta_indices[-eta_indices]
  }else{
    l_tot=p+(control$number_outcome_categories-1)
    beta_length=p
    eta_length=control$number_outcome_categories-1
    beta0_record = rep(0,beta_length) #Shared beta coefficients by site
    eta0_record = rep(0,eta_length) #Site specific eta coefficients by site
    eta_indices=(1:(control$number_outcome_categories-1))+p
    beta_indices=seq(1:l_tot)
    beta_indices=beta_indices[-eta_indices]
  }
  
  S_site_mat <- matrix(NA, K, beta_length)
  eta_mat=matrix(NA,K,eta_length)
  bbar_mat=matrix(NA, K, beta_length) #Will be the same for each site but just need to capture for this function
  nn=rep(0,K)
  for(i in 1:length(control$sites)){   
    site_i=control$sites[i]
    derivatives_i=pdaGet(paste0(site_i,'_derive'),config)
    S_site_mat[i,]=derivatives_i$S_site
    eta_mat[i,]=derivatives_i$eta
    nn[i]=derivatives_i$site_size
  }
  N=sum(nn)
  
  #Local Site
  local_site_name=config$site_id
  local_site=match(local_site_name,control$sites) #Local site index
  
  bbar=bbar_mat[local_site,]
  
  #Surrogate Likelihood Estimation
  x= as.matrix(ipdata[,2:ncol(ipdata)]) 
  y= as.matrix(ipdata[,1])
  
  if(control$ordinal_categories==FALSE){
    model="mlr"
  }else{
    model="polr"
  }
  
  #Summing up efficient score evaluated at each site for bbar as
  SN=colSums(S_site_mat)/N
  
  U_1=check_U(beta=bbar, betabar=bbar, eta_mat=eta_mat, site=local_site, x=x,y=y,k=local_site,eta_indices=eta_indices,model=model)
  
  S_es = function(beta){
    check_U(beta=beta, betabar=bbar, eta_mat=eta_mat, site=local_site, x=x,y=y,k=local_site,eta_indices=eta_indices,model=model)+ (SN - U_1) 
  }
  
  S_es_hessian=function(beta){
    numDeriv::jacobian(S_es, beta) 
  }
  
  S_es_euc=function(beta){
    sqrt(sum((check_U(beta=beta, betabar=bbar, eta_mat=eta_mat, site=local_site, x=x,y=y,k=local_site,eta_indices=eta_indices,model=model)+ (SN - U_1))^2)) 
  }
  
  root=optim(bbar,S_es_euc,method="BFGS", hessian=F, control = list(maxit=500,reltol=1e-8))
  
  ##Surrogate Efficient Score Results
  betanew = root$par
  
  Hessian=S_es_hessian(betanew)
  vcov=solve(-1*Hessian*N)
  beta_SE=sqrt(diag(vcov))
  
  
  if(control$ordinal_categories==FALSE){
    
    ses_res=list(btilde=betanew, 
                 btilde.se=beta_SE,
                 eta_mat=eta_mat,
                 eta_mat_theta=NULL,
                 site=config$site_id, 
                 site_size=n_site,
                 vcov=vcov)
  }else{
    
    eta_mat_alpha=matrix(NA,K,eta_length)
    for(i in 1:length(control$sites)){
      eta_theta=eta_mat[i,]
      eta_zeta=cumsum(c(eta_theta[1L], exp(eta_theta[-1L])))
      eta_mat_alpha[i,]=eta_zeta
    }
    
    
    ses_res=list(btilde=betanew, 
                 btilde.se=beta_SE,
                 eta_mat=eta_mat_alpha,
                 eta_mat_theta=eta_mat,
                 site=config$site_id, 
                 site_size=n_site,
                 vcov=vcov)
  }
  
  return(ses_res)
}


#' @useDynLib pda
#' @title PDA synthesize surrogate estimates from all sites, optional
#' 
#' @usage ODACATH.synthesize(ipdata,control,config)
#' @param ipdata local data in data frame
#' @param control pda control
#' @param config pda cloud configuration
#' 
#' @details Optional step-4: synthesize all the surrogate est btilde from each site, if step-3 from all sites is broadcasted
#'
#' @return  list(btilde=btilde, Vtilde=Vtilde)
#' @keywords internal
ODACATH.synthesize <- function(ipdata,control,config) {
  n=nrow(ipdata)  
  p=ncol(ipdata[,2:ncol(ipdata)])
  if(control$ordinal_categories==FALSE){
    px=(p)*(control$number_outcome_categories-1)
  }else{
    px=p
  }
  
  K <- length(control$sites)
  btilde_matrix=matrix(0,px,px)
  vcov_sum <- rep(0, px)  
  btilde.var_matrix=matrix(0,px,px)
  
  
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
