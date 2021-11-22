# Copyright 2021 Penn Computing Inference Learning (PennCIL) lab
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
DLM.steps <- c('initialize','estimate')
DLM.family <- 'gaussian'
# require(minqa)

# Distributed Linear Model
# 3 models: 
# - Linear model
# - Linear model with fixed site-effect
# - Linear model with random site-effect (Linear mixed model)

#' @useDynLib pda
#' @title DLM initialize
#' 
#' @usage DLM.initialize(ipdata,control,config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#' 
#' @references Yixin Chen, et al. (2006) Regression cubes with lossless compression and aggregation. 
#'    IEEE Transactions on Knowledge and Data Engineering, 18(12), pp.1585-1599. \cr
#'    (DLMM) Chongliang Luo, et al. (2020) Lossless Distributed Linear Mixed Model with Application to Integration of Heterogeneous Healthcare Data.  
#'    medRxiv, \doi{10.1101/2020.11.16.20230730}. \cr
#' @return init
#' @keywords internal
DLM.initialize <- function(ipdata,control,config){
  y = ipdata$outcome
  X = as.matrix(ipdata[,-'outcome']) # the first col is intercept?
  
  init = list(SiX = t(X) %*% X,
              SiXY = t(X) %*% y,
              SiY = sum(y^2),
              ni = length(y))
  return(init)
}

 


#' @useDynLib pda
#' @title PDA DLM estimation
#' 
#' @usage DLM.estimate(ipdata=NULL,control,config)
#' @param ipdata no need
#' @param control PDA control
#' @param config cloud configuration
#' 
#' @details DLM estimation: 
#' (1) Linear model, 
#' (2) Linear model with fixed effects, 
#' (3) Linear model with random effects (Linear mixed model)
#' @return  list(bhat, sebhat, sigmahat, uhat, seuhat) 
#' @keywords internal
DLM.estimate <- function(ipdata=NULL,control,config) {
  # data sanity check ... 
  
  K <- length(control$sites)
  SiXYZ <- list()
  for(site_i in control$sites){
    init_i <- pdaGet(paste0(site_i,'_initialize'),config)
    if(control$heterogeneity==T){
      # if vars with site-specific effects (fixed or random) are provided
      idx <- match(control$risk_factor_heterogeneity, control$risk_factor)
      SiXYZ[[site_i]] <- list(SiX  = init_i$SiX, SiXZ = as.matrix(init_i$SiX[,idx]), SiXY = init_i$SiXY,
                              SiZ  = as.matrix(init_i$SiX[idx,idx]), SiZY = init_i$SiXY[idx], 
                              SiY  = init_i$SiY, ni = init_i$ni)
    }else{
      SiXYZ[[site_i]] <- list(SiX  = init_i$SiX, SiXY = init_i$SiXY,
                              SiY  = init_i$SiY, ni = init_i$ni)
    }
  }
  
  px <- nrow(SiXYZ[[1]]$SiX)
  SX <- Reduce('+', lapply(SiXYZ, function(a) a$SiX))
  SXY <- Reduce('+', lapply(SiXYZ, function(a) a$SiXY))
  SY <- Reduce('+', lapply(SiXYZ, function(a) a$SiY))
  n <- Reduce('+', lapply(SiXYZ, function(a) a$ni))
  
  if(control$heterogeneity==F){
    ############## (1) Linear model ######################################
    # OLS
    bhat <- c(solve(SX, SXY))
    SSE <- SY - 2*sum(SXY*bhat) + t(bhat)%*%SX%*%bhat
    sigmahat <- c(sqrt(SSE/(n-px)))
    sebhat <- sqrt(diag(solve(SX))) * sigmahat
    
    names(sebhat) <- names(bhat) <- control$risk_factor
     
    res <- list(risk_factor=control$risk_factor, bhat=bhat, sebhat=sebhat, sigmahat=sigmahat)
  }else if(control$heterogeneity_effect=='fixed'){
    ############## (2) Linear model with fixed effects ##############
    # calculate total no. parameters p, warning if too many  
    q <- length(control$risk_factor_heterogeneity)
    p <- px + q*(K-1)  # lead site as the reference, so only K-1 levels of site effects
    if(p > 100) warning('Too many parameters (variables + variables_heterogeneity * sites > 100)!')
    ss <- setdiff(control$sites, control$lead_site)
    
    # block matrix
    BXZ <- matrix(0, p, p)
    BXY <- rep(0, p)
    BXZ[1:px, 1:px] <- SX
    BXY[1:px] <- SXY 
    for(si in seq_along(ss)){
      site_i <- ss[si]
      # indices for the block of ZtZ
      ib <- (px+(si-1)*q+1):(px+si*q)
      BXZ[1:px, ib] <- SiXYZ[[site_i]]$SiXZ
      BXZ[ib, 1:px] <- t(SiXYZ[[site_i]]$SiXZ)
      BXZ[ib, ib]   <- t(SiXYZ[[site_i]]$SiZ) 
      BXY[ib] <- SiXYZ[[site_i]]$SiZY 
    }
    
    # OLS
    bg = solve(BXZ, BXY)
    SSE <- SY - 2*sum(BXY*bg) + t(bg)%*%BXZ%*%bg
    sigmahat <- c(sqrt(SSE/(n-p)))
    sebg <- sqrt(diag(solve(BXZ))) * sigmahat
    
    bhat <- c(bg[1:px])
    sebhat <- sebg[1:px]
    uhat <- matrix(bg[-c(1:px)], nrow=K-1, byrow = T)
    seuhat <- matrix(sebg[-c(1:px)], nrow=K-1, byrow = T)
    
    names(sebhat) <- names(bhat) <- control$risk_factor
    row.names(uhat) <- row.names(seuhat) <- ss
    colnames(uhat) <- colnames(seuhat) <- control$risk_factor_heterogeneity
    
    res <- list(risk_factor=control$risk_factor, risk_factor_heterogeneity=control$risk_factor_heterogeneity, 
                bhat=bhat, sebhat=sebhat, sigmahat=sigmahat, uhat=uhat, seuhat=seuhat)
    # res <- lapply(res, function(a) as.data.frame(a))
  }else if(control$heterogeneity_effect=='random'){  # DLMM
    ############## (3) Linear model with random site effects (DLMM)  ##### 
    # LMM
    fit1 <- lmm.fit(SiXYZ = SiXYZ, pooled=F, reml=T, hessian=T) 
    # fit2 <- lmm.fit(Y=ipdata$outcome, X=as.matrix(ipdata[,-1]), Z=as.matrix(ipdata$`(Intercept)`),
    #                 id.site=LOS$site, pooled=T, reml=T, hessian=T) 
    
    bhat <- c(fit1$b)
    sebhat <- fit1$b.sd               # sd of fixed effect est
    uhat <- as.matrix(sapply(fit1$ui, function(a) a) )          # BLUP of random effects
    seuhat <- as.matrix(sqrt(sapply(fit1$varui_post, diag)))   # se of BLUPs
    sigmahat <- sqrt(fit1$s2)         # se of common error
    Vhat <- fit1$V                    # variance components: Var(ui)

    names(sebhat) <- names(bhat) <- control$risk_factor
    row.names(uhat) <- row.names(seuhat) <- control$sites
    colnames(uhat) <- colnames(seuhat) <- control$risk_factor_heterogeneity
    colnames(Vhat) <- row.names(Vhat) <- control$risk_factor_heterogeneity
   
    res <- list(risk_factor=control$risk_factor, risk_factor_heterogeneity=control$risk_factor_heterogeneity, 
                bhat=bhat, sebhat=sebhat, sigmahat=sigmahat, uhat=uhat, seuhat=seuhat, Vhat=Vhat)
    # res <- lapply(res, function(a) as.data.frame(a))
    # RJSONIO::fromJSON(RJSONIO::toJSON(res))
    # jsonlite::
    # rjson::    
  }
  
  return(res)
}

