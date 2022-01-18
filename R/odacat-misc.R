# Copyright 2022 Penn Computing Inference Learning (PennCIL) lab
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

## for internal use only: functions for DLMM and dPQL
# model.logL
# model.logL.d
# model.fit
# repar
# odacat.fit
# meta.fit

# require(MASS)
# require(numDeriv)

## (normalized neg) log-lik of proportional odds linear reg (polr) or 
## multinomial logistic reg (mlr) for ordinal outcome
#' @keywords internal
model.logL <- function(g, x, y , model='polr'){ # 'mlr'
  #POLR model is reparameterized to beta and theta
  # ordinal need to be coded as 1,2,3,...
  
  # if(missing(wt)) wt <- rep(1, length(y))
  # if(missing(offset)) offset <- rep(0, length(y))
  x <- as.matrix(x)
  n <- nrow(x)
  
  # for mlr, x need to contain intercept column!
  if(model=='mlr') x = cbind(1,x)
  
  px <- ncol(x)
  ind_px=1:px
  # lev <- levels(y)
  # if(length(lev) <= 2L) stop("response must have 3 or more levels")
  # y <- unclass(y)
  # q <- length(lev) - 1L
  # q = length(g)-px 
  q <- length(unique(y)) - 1
  ind_q <- seq_len(q) 
  # expit
  pfun <- function(a) 1/(1+exp(-a))
  
  if(q <= 1L) stop("response must have 3 or more levels")
  
  if(model=='polr'){
    theta <- g[px + ind_q]
    gamm <- c(-Inf, cumsum(c(theta[1L], exp(theta[-1L]))), Inf)
    # eta <- offset
    # if(px ) eta <- eta + drop(x %*% g[ind_px])
    if(px ) eta <- drop(x %*% g[ind_px])
    z1 <- pmin(100, gamm[y+1L] - eta)
    z2 <- pmax(-100, gamm[y] - eta)
    pr <- pfun(z1) - pfun(z2)
    # logL <- if (all(pr > 0)) -sum(wt * log(pr)) else Inf
    logL <- if (all(pr > 0)) -sum(log(pr)) else Inf
  }else if(model=='mlr'){  
    bm = cbind(matrix(g, px, q),0)
    Xb = x %*% bm
    li=sapply(1:nrow(x),function(i){Xb[i,y[i]]})-log(apply(exp(Xb),1,sum))
    #li = diag(Xb[,y]) - log(apply(exp(Xb),1,sum)) slow for large datasets
    logL =  - sum(li)
  }
  return(logL/n)
}

## gradient and hessian of polr.logL
#' @keywords internal
model.logL.d <- function(g, x, y, model='polr', grad=T, hess=T){
  # if(missing(wt)) wt <- rep(1, length(y))
  # if(missing(offset)) offset <- rep(0, length(y))
  ll <- function(g) model.logL(g, x=x, y=y, model=model)
  logL <- ll(g)
  logL.grad = logL.hess = NULL
  
  if(grad){
    logL.grad <- numDeriv::grad(ll, g)  
    
    # using numDeriv is convenient, but much slower (~20 times) than manual derivation
    ### manually derive the gradient...    
    #   jacobian <- function(theta) { ## dgamma by dtheta matrix
    #     k <- length(theta)
    #     etheta <- exp(theta)
    #     mat <- matrix(0 , k, k)
    #     mat[, 1L] <- rep(1, k)
    #     for (i in 2L:k) mat[i:k, i] <- etheta[i]
    #     mat
    #   }
    #   
    #   p1 <- dfun(z1); p2 <- dfun(z2)
    #   g1 <- if(px ) t(x) %*% (wt*(p1 - p2)/pr) else numeric()
    #   xx <- .polrY1*p1 - .polrY2*p2
    #   g2 <- - t(xx) %*% (wt/pr)
    #   g2 <- t(g2) %*% jacobian(theta)
    #   logL.grad = if(all(pr > 0)) c(g1, g2) else rep(NA_real_, px+q)
  }
  
  if(hess){ 
    logL.hess <- numDeriv::hessian(ll, g)
  }
  
  return(list(logL=logL, logL.grad=logL.grad, logL.hess=logL.hess))
}

## fit polr / mlr
#' @keywords internal
model.fit <- function(x, y, model='polr', b0){
  # if(missing(wt)) wt <- rep(1, length(y))
  # if(missing(offset)) offset <- rep(0, length(y))
  n <- length(y)
  x <- as.matrix(x)
  px <- ncol(x)
  q <- length(unique(y)) - 1
  if(q <= 1L) stop("response must have 3 or more levels")
  # ind_q <- seq_len(q)
  pp = ifelse(model=='polr', px+q, (px+1)*q)
  if(missing(b0)) b0 = rep(0, pp) 
  
  # ll <- function(g) polr.logL(g, x, y, wt, offset, method) 
  res <- optim(b0, model.logL, x=x, y=y, model=model, method="BFGS", hessian=T, control = list(maxit=1000,reltol=1e-6))

  tt = list()
  
  if(model=='polr'){ #Returns back to zeta scale for intercept
    tt <- repar(g=res$par, g.covar=solve(res$hessian*n), px=px)
    if(px) names(tt$beta) <- colnames(x)
  }
  
  # names(tt$zeta) <- paste(lev[-length(lev)], lev[-1L], sep="|")
  # names(tt$theta) <- paste(lev[-length(lev)], lev[-1L], sep="|")
  
  tt$res <- res
  tt$deviance <- 2 * n * res$value
  tt$model <- model
  if(model=='mlr'){
    tt$var_cov=solve(res$hessian*n)
  }
  return(tt) 
}

## reparar from gamma to beta, zeta, and calculate var(zeta)
#' @keywords internal
repar <- function(g, g.covar, px){  
  q <- length(g) - px
  ind_q <- seq_len(q)
  beta = g[1:px]
  theta <- g[px + ind_q]
  zeta <- cumsum(c(theta[1L], exp(theta[-1L])))
  beta.var <- theta.var <- zeta.var <- rep(NA, q)
  if(!missing(g.covar)){
    H <- g.covar  # solve(g.hess * n)
    beta.var <- diag(H)[seq_len(px)]
    theta.var <- diag(H)[px+seq_len(q)]
    zeta.var[1] <- diag(H)[px+1]
    for(iq in 2:q){                           # calculate var(zeta[iq]) using delta method
      Hi = H[px+seq_len(iq), px+seq_len(iq)]  # diag(H)[px+seq_len(iq)]
      di = c(1, exp(theta[2:iq]))             # d(theta_1 + e^theta_2 + ...+e^theta_iq) / d(theta)
      zeta.var[iq] <- t(di) %*% Hi %*% di
    }
  }
  
  return(list(beta = beta, theta = theta, zeta = zeta, 
              beta.var = beta.var, theta.var=theta.var, zeta.var = zeta.var))
}

#Beta bar in terms of beta and theta for model
#' @keywords internal
odacat.fit <- function(x, y, 
                       # mydata, # col1=site, col2=y, col3+=X 
                       # meta_est_list,
                      id.local = NULL, 
                      # method='logistic',
                      # init_est = 'meta', 
                      beta_bar,
                      logL_N_D1_beta_bar,
                      logL_N_D2_beta_bar,
                      nn,
                      model='polr',
                      verbose = F,  
                      optim_control = list(maxit=10000)){
  
  
  # ODAM1 and ODAM2
  sum_D1_beta_bar <- colSums(logL_N_D1_beta_bar * nn) / sum(nn)
  sum_D2_beta_bar <- apply(logL_N_D2_beta_bar * nn, c(2,3), sum) / sum(nn)
  
  # di <- as.matrix(mydata)
  grad_LN_L1_beta_bar <- sum_D1_beta_bar - logL_N_D1_beta_bar[id.local,]
  hess_LN_L1_beta_bar <- sum_D2_beta_bar - logL_N_D2_beta_bar[id.local,,]
  
  # surrogate log-L 
  Ltilde1 <- function(g) model.logL(g, x=x, y=y, model=model) + sum(g * grad_LN_L1_beta_bar)
  Ltilde2 <- function(g) Ltilde1(g) + 1/2 * t(c(g-beta_bar)) %*% hess_LN_L1_beta_bar %*% c(g-beta_bar)
  
  sol.odacat2 <- optim(par = beta_bar, fn = Ltilde2, hessian = F, control = optim_control)
  sol.odacat2.hess = numDeriv::hessian(Ltilde2, sol.odacat2$par)
  
  
  return(list(beta_bar = beta_bar, 
              sol.odacat2=sol.odacat2,
              sol.odacat2.hess=sol.odacat2.hess))
  
}

#meta.fit provides parameters in terms of beta and theta (reparameterized intercept)
#' @keywords internal
meta.fit <- function(p, # ipdata,
                     control,config){
  # n=nrow(ipdata)
  # p=ncol(ipdata[,2:ncol(ipdata)])
  if(control$ordinal_categories==FALSE){ #MLR
  px=(p+1)*(control$number_outcome_categories-1)
  }else{ #POLR
  px=p+(control$number_outcome_categories-1) 
  }
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
  
 return(list(bbar=betameta,vmeta=vmeta))
}


