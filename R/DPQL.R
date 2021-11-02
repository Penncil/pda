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
# DPQL.steps <- c('initialize', 'derive', 'estimate')
# DPQL.family <- 'binomial' # family as in glm
# require(minqa)

# distributed Penalized Quasi-Likelihood (dPQL) algorithm  
# for fitting multi-center Generalized Linear Mixed Model (GLMM)

#' @useDynLib pda
#' @title DPQL initialize
#' 
#' @usage DPQL.initialize(ipdata,control,config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#' @details To initialize, fit glm at each individual site and send the estimated effect size and variances to the lead site. 
#'          This step may be optional if we just use zero's as initial effect sizes to start the PQL algorithm.
#' @return init
#' @keywords internal
DPQL.initialize <- function(ipdata,control,config){ 
  fit_i <- glm(outcome ~ ., data=ipdata, family=control$family)
  
  init <- list(bhat_i = fit_i$coef,
               Vhat_i = summary(fit_i)$coef[,2]^2,   
               site = config$site_id,
               site_size = nrow(ipdata))
  return(init)
}


#' @useDynLib pda
#' @title DPQL derive
#' 
#' @usage DPQL.derive(ipdata,control,config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#' @details This step calculated the intermediate aggregated data (XtWX, XtWY, and YtWY) for each site. 
#'           May need to be iterated several times until prespecified rounds are met.
#' @references Chongliang Luo, et al. (2021) dPQL: a lossless distributed algorithm for generalized linear mixed model 
#'                   with application to privacy-preserving hospital profiling. medRxiv, \doi{10.1101/2021.05.03.21256561}. \cr
#'            Chongliang Luo, et al. (2020) Lossless Distributed Linear Mixed Model with Application to Integration of Heterogeneous Healthcare Data.  
#'                  medRxiv, \doi{10.1101/2020.11.16.20230730}. \cr
#' @return list(SiX, SiXY, SiY, ni)
#' @keywords internal
DPQL.derive <- function(ipdata,control,config){
  ## distributed PQL (DPQL)  
  ## collaborative site calculate aggregated data 
  ## t(Xi) %*% Wi %*% Xi, t(Xi) %*% Wi %*% Yi, t(Yi) %*% Wi %*% Yi
  family <- control$family 
  if (is.character(family)) 
    family <- get(family)
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  Y <- ipdata$outcome
  X <- as.matrix(ipdata[,-'outcome'])  
  idx <- match(control$risk_factor_heterogeneity, control$risk_factor)
  Z <- as.matrix(X[,idx])
  # if(is.null(weights)) weights <- rep(1, length(Y)) 
  # if(is.null(offset)) offset <- rep(0, length(Y))
  
  ## intermediate fixed and random effects
  # qz <- ncol(Z); px <- ncol(X)
  # if(is.null(ranef.init)) ranef.init=rep(0,qz)
  # if(is.null(fixef.init)) fixef.init=rep(0,px)
  if(is.null(control$bhat)) fixef.init <- rep(0,ncol(X)) else fixef.init <- control$bhat
  if(is.null(control$uhat)) ranef.init <- rep(0,ncol(Z)) else ranef.init <- control$uhat[control$sites==config$site_id,] 
 
  # # get b_meta as initial bbar
  # bhat <- rep(0, px)
  # vbhat <- rep(0, px)     # cov matrix?
  # for(site_i in control$sites){
  #   init_i <- pdaGet(paste0(site_i,'_initialize'),config)
  #   bhat = rbind(bhat, init_i$bhat_i)
  #   vbhat = rbind(vbhat, init_i$Vhat_i)
  # }
  # bhat = bhat[-1,]
  # vbhat = vbhat[-1,]
  # 
  # #estimate from meta-analysis
  # betameta = apply(bhat/vbhat,2,function(x){sum(x, na.rm = T)})/apply(1/vbhat,2,function(x){sum(x, na.rm = T)})
  # vmeta = 1/apply(1/vbhat,2,function(x){sum(x, na.rm = T)})
 
  ## weight and pseudo outcome 
  # w <- weights
  w <- rep(1, length(Y)) 
  eta <- c(X%*%fixef.init + Z %*% ranef.init) 
  mu <- family$linkinv(eta)
  mu.eta.val <- family$mu.eta(eta)
  wz <- w * mu.eta.val^2 / family$variance(mu) 
  zz <- eta + (Y-mu)/mu.eta.val # - offset 
  
  ## aggregated data: (XtWX, XtWY, and YtWY)
  SiXYZ <- list(SiX  = t(X*wz) %*% X, 
                SiXY = t(X*wz) %*% zz, 
                SiY  = sum(zz^2*wz),
                ni = length(Y)) 
  
  return(SiXYZ)
}


#' @useDynLib pda
#' @title PDA DPQL estimation
#' 
#' @usage DPQL.estimate(ipdata=NULL,control,config)
#' @param ipdata no need
#' @param control PDA control
#' @param config cloud configuration
#' 
#' @details DPQL estimation: (iterative) weighted DLMM using AD from all sites
#' @references Chongliang Luo, et al. (2021) dPQL: a lossless distributed algorithm for generalized linear mixed model 
#'                   with application to privacy-preserving hospital profiling. medRxiv, \doi{10.1101/2021.05.03.21256561}. \cr
#'            Chongliang Luo, et al. (2020) Lossless Distributed Linear Mixed Model with Application to Integration of Heterogeneous Healthcare Data.  
#'                 medRxiv, \doi{10.1101/2020.11.16.20230730}. \cr
#' @return  list(risk_factor, risk_factor_heterogeneity, bhat, sebhat, uhat, seuhat, Vhat) 
#' @keywords internal
DPQL.estimate <- function(ipdata=NULL,control,config) {
  # data sanity check ... 
  
  K <- length(control$sites)
  # iround <- as.numeric(gsub('[^[:digit:]]', '', control$step) )
  SiXYZ <- list()
  for(site_i in control$sites){
    init_i <- pdaGet(paste0(site_i,'_derive_', control$round), config)
    idx <- match(control$risk_factor_heterogeneity, control$risk_factor)
    SiXYZ[[site_i]] <- list(SiX  = init_i$SiX, SiXZ = as.matrix(init_i$SiX[,idx]), SiXY = init_i$SiXY,
                            SiZ  = as.matrix(init_i$SiX[idx,idx]), SiZY = init_i$SiXY[idx], 
                            SiY  = init_i$SiY, ni = init_i$ni)
  }

  ## weighted DLMM 
  fit1 <- lmm.fit(SiXYZ = SiXYZ, pooled=F, reml=T, hessian=T)  
  
  bhat <- c(fit1$b)
  sebhat <- fit1$b.sd               # sd of fixed effect est
  uhat <- as.matrix(sapply(fit1$ui, function(a) a))          # BLUP of random effects
  seuhat <- as.matrix(sqrt(sapply(fit1$varui_post, diag)))   # se of BLUPs
  # sigmahat <- sqrt(fit1$s2)         # se of common error
  Vhat <- fit1$V                    # variance components: Var(ui)
  
  names(sebhat) <- names(bhat) <- control$risk_factor
  row.names(uhat) <- row.names(seuhat) <- control$sites
  colnames(uhat) <- colnames(seuhat) <- control$risk_factor_heterogeneity
  colnames(Vhat) <- row.names(Vhat) <- control$risk_factor_heterogeneity
  
  res <- list(risk_factor=control$risk_factor, 
              risk_factor_heterogeneity=control$risk_factor_heterogeneity, 
              bhat=bhat, sebhat=sebhat, #sigmahat=sigmahat, 
              uhat=uhat, seuhat=seuhat, Vhat=Vhat)
    
  return(res)
}

