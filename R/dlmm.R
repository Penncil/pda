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

## for internal use only: functions for DLMM and dPQL
# lmm.profile
# lmm.profile0
# lmm.fit
# lmm.lrt
# lmm.get.summary
# lmm.get.summary.cut
# lmm.fit.partial
# dlmm
# glmmDPQL.fit

 
# require(Rcpp)
# Rcpp::sourceCpp('dlmm.cpp')
# require(minqa)

#' @keywords internal
lmm.profile <- function(V, s2,             # var components para
                        pooled = F, reml = T, 
                        Y, X, Z, id.site, weights=NULL,  # if pooled ipd is available
                        SiXYZ,             # if pooled ipd is not available, use summary stats
                        rcpp = F){
  ## lmm profile likelihood w.r.t. variance components
  if(pooled == T){
    id.site.uniq <- unique(id.site)
    px <- ncol(X)
    pz <- ncol(Z)
  }else{ 
    id.site.uniq <- names(SiXYZ)
    px <- ncol(SiXYZ[[1]]$SiX)
    pz <- ncol(SiXYZ[[1]]$SiXZ)
  }
  
  # allows assuming the same residual var across sites
  if(length(s2) == 1) {
    s2 <- rep(s2, length(id.site.uniq))
  }
  
  if(rcpp == F){    
    lpterm1 = lpterm2 = remlterm = 0
    remlterm.i = rep(NA, length(id.site.uniq))
    
    bterm1 <- matrix(0, px, px)
    bterm2 <- rep(0, px)
    Vinv <- solve(V)
    Wi <- ui <- varui <- varui_post <- list()  # save this for each subject for BLUP
    
    for(ii in seq_along(id.site.uniq)){
      si <- id.site.uniq[ii]
      if(pooled == T){
        SiXYZ <- lmm.get.summary(Y, X, Z, weights, id.site)
      } #else{
      SiX  <- SiXYZ[[si]]$SiX
      SiXZ <- SiXYZ[[si]]$SiXZ
      SiXY <- SiXYZ[[si]]$SiXY
      SiZ <- SiXYZ[[si]]$SiZ
      SiZY <- SiXYZ[[si]]$SiZY
      SiY  <- SiXYZ[[si]]$SiY
      ni <- SiXYZ[[si]]$ni
      # }
      
      s2i <- s2[ii]  # sigma_i^2
      
      tmp <- log(det(diag(1, pz) + SiZ %*% V / s2i))
      
      if(is.na(tmp)) cat(diag(V), '...', s2i, '...', V[1,], '\n')
      # logdet <- ni * log(s2i) + log(det(diag(1, pz) + SiZ %*% V / s2i))
      logdet <- ni * log(s2i) + log(det(diag(1, pz) + SiZ %*% V / s2i))
      
      # log(max(1e-14, det(diag(1,pz)+SiZ%*%V/s2i)))
      lpterm1 <- lpterm1 + logdet
      
      Wi[[ii]] <- solve(s2i * Vinv + SiZ)
      bterm1 <- bterm1 + (SiX - SiXZ %*% Wi[[ii]] %*% t(SiXZ)) / s2i
      bterm2 <- bterm2 + (SiXY - SiXZ %*% Wi[[ii]] %*% SiZY) / s2i
      lpterm2 <- lpterm2 + (SiY - t(SiZY) %*% Wi[[ii]] %*% SiZY) / s2i
      
      if(reml == T){
        # tmp = log(det((SiX-SiXZ%*%Wi[[ii]]%*%t(SiXZ))/s2i))
        # if(is.na(tmp)) cat(Wi[[ii]], '...', s2i, '\n')
        remlterm.i[ii] <- log(abs(det((SiX - SiXZ %*% Wi[[ii]] %*% t(SiXZ))/s2i))) # abs() to make sure positivity
        # remlterm <- remlterm + log(det((SiX - SiXZ %*% Wi[[ii]] %*% t(SiXZ)) / s2i))
      }
    }
    
    b <- solve(bterm1, bterm2)
    if(reml == T){
      remlterm = sum(remlterm.i[is.finite(remlterm.i)])
      lp <- - (lpterm1 + lpterm2 - 2 * sum(bterm2 * b) + t(b) %*% bterm1 %*% b + remlterm) / 2
    }else{
      lp <- - (lpterm1 + lpterm2 - 2 * sum(bterm2 * b) + t(b) %*% bterm1 %*% b) / 2
    }
    
    lk <- -(lpterm1 + lpterm2 - 2 * sum(bterm2 * b) + t(b) %*% bterm1 %*% b) / 2 # To be used in calculating mAIC
    
    # Loop to estimate ui
    # ui = V Z_i^{T} Sigma_{i}^{-1} (Y_i - X_i \beta)
    for(ii in seq_along(id.site.uniq)){ # re-call this loop
      si <- id.site.uniq[ii]
      if(pooled == T){
        Xi <- X[id.site == si, ]
        Zi <- Z[id.site == si, ]
        Yi <- Y[id.site == si]
        SiX  <- t(Xi) %*% Xi
        SiXZ <- t(Xi) %*% Zi
        SiXY <- t(Xi) %*% Yi
        SiZ  <- t(Zi) %*% Zi
        SiZY <- t(Zi) %*% Yi
        SiY  <- sum(Yi ^ 2)
        ni <- sum(id.site == si)
      }else{
        SiX  <- SiXYZ[[si]]$SiX
        SiXZ <- SiXYZ[[si]]$SiXZ
        SiXY <- SiXYZ[[si]]$SiXY
        SiZ <- SiXYZ[[si]]$SiZ
        SiZY <- SiXYZ[[si]]$SiZY
        SiY <- SiXYZ[[si]]$SiY
        ni <- SiXYZ[[si]]$ni
      }
      
      s2i <- s2[ii]  # sigma_i^2
      uiterm1 <- (SiZY - SiZ %*% Wi[[ii]] %*% SiZY) / s2i
      uiterm2 <-  ((t(SiXZ) - SiZ %*% Wi[[ii]] %*% t(SiXZ)) / s2i) %*% as.matrix(b) 
      ui[[ii]] <- V %*% as.numeric(uiterm1 - uiterm2)
      
      vterm1 <- V %*% ((SiZ - SiZ %*% Wi[[ii]] %*% SiZ) / s2i) %*% V
      vterm2 <- V %*% ((t(SiXZ) - SiZ %*% Wi[[ii]] %*% t(SiXZ)) / s2i)
      # varui[[ii]] <- V - vterm1 + (vterm2 %*% t(vterm2) / lpterm1)
      # varui_post[[ii]] <- vterm1 - (vterm2 %*% t(vterm2) / lpterm1)
      # 20210112:
      varui[[ii]] <- V - vterm1 + (vterm2 %*% solve(bterm1, t(vterm2)))
      varui_post[[ii]] <- vterm1 - (vterm2 %*% solve(bterm1, t(vterm2)))
    }
    
    
    res <- list(lp = lp, b = b, ui = ui, lk = lk, varui = varui, varui_post = varui_post,
                allterms = list(lpterm1 = lpterm1,
                                lpterm2 = lpterm2,
                                remlterm = remlterm,
                                remlterm.i = remlterm.i,
                                bterm1 = bterm1,
                                bterm2 = bterm2))
  } else{ ## CHECK THIS 
    res <- NULL # LMM_Profile(SiXYZ = SiXYZ, V = V, s2 = s2, reml = reml)  
  }
  # SiY - t(SiZY)%*%Wi%*%SiZY - 2*(t(SiXY)-t(SiZY)%*%Wi%*%t(SiXZ))%*%b + t(b)%*%(SiX-SiXZ%*%Wi%*%t(SiXZ))%*%b
  return(res) 
}


#' @keywords internal
lmm.profile0 <- function(V,                # var components para, = original V / s2    # s2, 
                         pooled = F, reml = T, 
                         Y, X, Z, id.site, weights=NULL,  # if pooled ipd is available
                         SiXYZ,             # if pooled ipd is not available, use summary stats
                         rcpp = F){
  ## further profile out the residual var s2, used if common.s2=T (deemed as usual situation)
  if(pooled == T){
    id.site.uniq <- unique(id.site)
    px <- ncol(X)
    pz <- ncol(Z)
  }else{ 
    id.site.uniq <- names(SiXYZ)
    px <- ncol(SiXYZ[[1]]$SiX)
    pz <- ncol(SiXYZ[[1]]$SiXZ)
  }
  
  if(rcpp == F){    
    lpterm1 = lpterm2 = remlterm = 0
    bterm1 <- matrix(0, px, px)  # sum_i Xi' \Sigma_i^-1 Xi
    bterm2 <- rep(0, px)         # sum_i Xi' \Sigma_i^-1 Yi
    Vinv <- solve(V)
    Wi <- ui <- varui <- varui_post <- list()  # save this for each subject for BLUP
    
    N <- 0
    for(ii in seq_along(id.site.uniq)){
      si <- id.site.uniq[ii]
      if(pooled == T){
        SiXYZ <- lmm.get.summary(Y, X, Z, weights, id.site)
      } # else{
      SiX  <- SiXYZ[[si]]$SiX
      SiXZ <- SiXYZ[[si]]$SiXZ
      SiXY <- SiXYZ[[si]]$SiXY
      SiZ <- SiXYZ[[si]]$SiZ
      SiZY <- SiXYZ[[si]]$SiZY
      SiY  <- SiXYZ[[si]]$SiY
      ni <- SiXYZ[[si]]$ni
      # }
      N <- N + ni
      
      # tmp <- log(det(diag(1, pz) + SiZ %*% V))  # improve
      logdet <- log(det(diag(1, pz) + SiZ %*% V))   # -log(det(diag(wti))) omitted as wti are the fixed weights
      lpterm1 <- lpterm1 + logdet
      
      Wi[[ii]] <- solve(Vinv + SiZ)
      bterm1 <- bterm1 + (SiX - SiXZ %*% Wi[[ii]] %*% t(SiXZ)) # / s2i
      bterm2 <- bterm2 + (SiXY - SiXZ %*% Wi[[ii]] %*% SiZY) # / s2i
      lpterm2 <- lpterm2 + (SiY - t(SiZY) %*% Wi[[ii]] %*% SiZY) #/ s2i
    }
    
    b <- solve(bterm1, bterm2)
    # quadratic term: = (Yi-Xi*b)'\Sigma_i^-1 (Yi-Xi*b)
    qterm <- as.numeric(lpterm2 - 2 * sum(bterm2 * b) + t(b) %*% bterm1 %*% b ) 
    if(reml == T){
      remlterm = log(det(bterm1))
      s2 <- qterm / (N-px)
      # lp <- - (lpterm1 + lpterm2 - 2 * sum(bterm2 * b) + t(b) %*% bterm1 %*% b + remlterm) / 2
      lp <- - (lpterm1 + (1 + log(qterm * 2 * pi / (N - px))) * (N - px)) / 2     # Bates2015JSS eq(42)
    }else{
      s2 <- qterm / N
      lp <- - (lpterm1 + (1 + log(qterm * 2 * pi / N)) * N) / 2               # Bates2015JSS eq(35)
    }
    
    # lk <- -(lpterm1 + lpterm2 - 2 * sum(bterm2 * b) + t(b) %*% bterm1 %*% b) / 2 # To be used in calculating mAIC
    lk <- - (lpterm1 + (1+log(qterm*2*pi/N))*N) / 2 
    
    # Loop to estimate ui
    # ui = V Z_i^{T} Sigma_{i}^{-1} (Y_i - X_i \beta)
    for(ii in seq_along(id.site.uniq)){ # re-call this loop
      si <- id.site.uniq[ii]
      
      SiX  <- SiXYZ[[si]]$SiX
      SiXZ <- SiXYZ[[si]]$SiXZ
      SiXY <- SiXYZ[[si]]$SiXY
      SiZ <- SiXYZ[[si]]$SiZ
      SiZY <- SiXYZ[[si]]$SiZY
      SiY <- SiXYZ[[si]]$SiY
      ni <- SiXYZ[[si]]$ni 
      
      s2i <- s2   # [ii]  # sigma_i^2
      uiterm1 <- (SiZY - SiZ %*% Wi[[ii]] %*% SiZY) / s2i
      uiterm2 <-  ((t(SiXZ) - SiZ %*% Wi[[ii]] %*% t(SiXZ)) / s2i) %*% as.matrix(b) 
      ui[[ii]] <- V %*% as.numeric(uiterm1 - uiterm2) * s2i #
      
      vterm1 <- V %*% ((SiZ - SiZ %*% Wi[[ii]] %*% SiZ) / s2i) %*% V * (s2i ^ 2)    #
      vterm2 <- V %*% ((t(SiXZ) - SiZ %*% Wi[[ii]] %*% t(SiXZ)) / s2i)  * s2i   #
      varui[[ii]] <- (V * s2i) - vterm1 + (vterm2 %*% solve(bterm1, t(vterm2)))     # 20210111
      varui_post[[ii]] <- vterm1 - (vterm2 %*% solve(bterm1, t(vterm2)))            # 20210111
    }
    
    res <- list(lp = lp, b = b, 
                s2 = s2,
                ui = ui, lk = lk, varui = varui, varui_post = varui_post,  
                allterms = list(lpterm1 = lpterm1,
                                lpterm2 = lpterm2,
                                qterm = qterm,
                                remlterm = remlterm,
                                # remlterm.i = remlterm.i,
                                bterm1 = bterm1,
                                bterm2 = bterm2))
  } else{ ## CHECK THIS 
    res <- NULL # LMM_Profile(SiXYZ = SiXYZ, V = V, s2 = s2, reml = reml)  
  }
  return(res) 
}


#' @keywords internal
lmm.fit <- function(Y = NULL, X = NULL, Z = NULL, id.site = NULL, weights = NULL,
                    pooled = F, reml = T, 
                    common.s2 = T,      # common residual var across sites
                    SiXYZ = list(),
                    corstr = 'independence', # 'exchangeable', 'ar1', 'unstructured'),
                    mypar.init = NULL,
                    hessian = F,
                    verbose){
  ## fit lmm distributed
  if(pooled == T){
    id.site.uniq <- unique(id.site)
    px <- ncol(X)
    pz <- ncol(Z)
    K <- length(id.site.uniq)
    SiXYZ <- lmm.get.summary(Y, X, Z, weights, id.site)
  }else{ 
    id.site.uniq <- names(SiXYZ)
    px <- ncol(SiXYZ[[1]]$SiX)
    pz <- ncol(SiXYZ[[1]]$SiXZ)
    K <- length(SiXYZ)
  }
  
  ## 20200629: now further profile out s2 if common.s2=T
  if(common.s2 == T){
    ns <- 1 # number of s2 para
    fn <- function(mypar){ 
      if(corstr == 'independence'){
        V <- diag(mypar[1 : pz], pz)
        # V <- diag(exp(mypar[1 : pz]), pz)
        # s2 <- exp(mypar[-c(1 : pz)])
      }else if(corstr == 'exchangeable'){
        V = diag(sqrt(mypar[1 : pz])) %*% (matrix(mypar[pz + 1], pz, pz) + diag(1 - mypar[pz + 1], pz)) %*% diag(sqrt(mypar[1 : pz]))
        # s2 <- mypar[-c(1 : (pz + 1))]   
      }else if(corstr == 'unstructured'){
        V = matrix(0, pz, pz)
        diag(V) <- mypar[1 : pz]
        V[lower.tri(V)] <- V[upper.tri(V)] <- mypar[(pz + 1) : (pz * (pz + 1) / 2)]
        # s2 <- mypar[-c(1 : (pz*(pz+1)/2))]
      }
      return(-lmm.profile0(V, pooled=F, reml, Y, X, Z, id.site, weights, SiXYZ)$lp)
    }
    
    if(is.null(mypar.init)){
      if(corstr == 'independence'){
        mypar.init <- rep(0.5, pz)
        # mypar.init <- log(c(rep(0.5, pz)))
      }else if(corstr == 'exchangeable'){
        mypar.init <- c(rep(0.5, pz), 0.1 )
      }else if(corstr == 'unstructured'){
        mypar.init <- c(rep(0.5, pz), rep(0.1, pz * (pz - 1) / 2) )
      }
      cat('default mypar.init (var comp) = ', mypar.init, '\n')
    }
    
    # res <- optim(mypar.init, fn, hessian = hessian)
    res <- minqa::bobyqa(mypar.init, fn, lower=rep(1e-6, pz), control=list(maxfun=1e5))
    # res <- nloptwrap(mypar.init, fn, lower=rep(1e-6, pz), upper = rep(1e6, pz), control=list(maxfun=1e5))
    
    mypar <- res$par
    if(corstr == 'independence'){
      V <- diag(mypar[1 : pz], pz)
      # V <- diag(exp(mypar[1 : pz]), pz)
      # s2 <- exp(mypar[- c(1 : pz)])
    }else if(corstr == 'exchangeable'){
      V <- diag(sqrt(mypar[1 : pz])) %*% (matrix(mypar[pz + 1], pz, pz) + diag(1 - mypar[pz + 1], pz)) %*% diag(sqrt(mypar[1 : pz]))
      # s2 <- mypar[- c(1 : (pz + 1))]
    }else if(corstr == 'unstructured'){
      V <- matrix(0, pz, pz)
      diag(V) <- mypar[1 : pz]
      V[lower.tri(V)] <- V[upper.tri(V)] <- mypar[(pz + 1) : (pz * (pz + 1) / 2)]
      # s2 <- mypar[-c(1 : (pz * (pz + 1) / 2))]
      # error('corstr=="unstructured" not yet implemented')
    }
    
    res.profile <- lmm.profile0(V = V, pooled=F, reml, Y, X, Z, id.site, weights, SiXYZ)
    s2 <- res.profile$s2
    V <- V * s2             # scale back
  }else{  # if common.s2=F, can't profile out s2 vector
    ns <- K 
    fn <- function(mypar){
      if(corstr == 'independence'){
        V <- diag(mypar[1 : pz], pz)
        s2 <- exp(mypar[-c(1 : pz)])
      }else if(corstr == 'exchangeable'){
        V = diag(sqrt(mypar[1 : pz])) %*% (matrix(mypar[pz + 1], pz, pz) + diag(1 - mypar[pz + 1], pz)) %*% diag(sqrt(mypar[1 : pz]))
        s2 <- mypar[-c(1 : (pz + 1))]   
      }else if(corstr == 'unstructured'){
        V = matrix(0, pz, pz)
        diag(V) <- mypar[1 : pz]
        V[lower.tri(V)] <- V[upper.tri(V)] <- mypar[(pz + 1) : (pz * (pz + 1) / 2)]
        s2 <- mypar[-c(1 : (pz*(pz+1)/2))]
        # error('corstr=="unstructured" not yet implemented')
      }
      return(-lmm.profile(V, s2, pooled=F, reml, Y, X, Z, id.site, weights, SiXYZ)$lp)
    }
    
    if(is.null(mypar.init)){
      if(corstr == 'independence'){
        mypar.init <- c(rep(0.5, pz), rep(0.5, ns))
      }else if(corstr == 'exchangeable'){
        mypar.init <- c(rep(0.5, pz), 0.1, rep(0.5, ns))
      }else if(corstr == 'unstructured'){
        mypar.init <- c(rep(0.5, pz), rep(0.1, pz * (pz - 1) / 2), rep(0.5, ns))
      }
      cat('default mypar.init (var comp) = ', mypar.init, '\n')
    }
    
    # res <- optim(mypar.init, fn, hessian = hessian)
    res <- minqa::bobyqa(mypar.init, fn, lower=rep(1e-6, length(mypar.init)), control=list(maxfun=1e5))
    # res <- nloptwrap(mypar.init, fn, lower=rep(1e-6, length(mypar.init)), 
    #                  upper =rep(1e6, length(mypar.init)),  control=list(maxfun=1e5))
    
    mypar <- res$par
    if(corstr == 'independence'){
      V <- diag(mypar[1 : pz], pz)
      s2 <- mypar[- c(1 : pz)]
    }else if(corstr == 'exchangeable'){
      V <- diag(sqrt(mypar[1 : pz])) %*% (matrix(mypar[pz + 1], pz, pz) + diag(1 - mypar[pz + 1], pz)) %*% diag(sqrt(mypar[1 : pz]))
      s2 <- mypar[- c(1 : (pz + 1))]
    }else if(corstr == 'unstructured'){
      V <- matrix(0, pz, pz)
      diag(V) <- mypar[1 : pz]
      V[lower.tri(V)] <- V[upper.tri(V)] <- mypar[(pz + 1) : (pz * (pz + 1) / 2)]
      s2 <- mypar[-c(1 : (pz * (pz + 1) / 2))]
      # error('corstr=="unstructured" not yet implemented')
    }
    
    res.profile <- lmm.profile(V = V, s2 = s2, pooled, reml, Y, X, Z, id.site, SiXYZ) 
  }
  
  
  
  
  
  
  ## New added 
  ## Inference (Wald test statistic)
  vd <- diag(solve(res.profile$allterms$bterm1))
  if(common.s2==T) vd <- diag(solve(res.profile$allterms$bterm1 / s2))  # scale back
  wald <- res.profile$b / sqrt(vd)
  
  ## 95% CI for fixed effects  
  lb <- res.profile$b -  1.96 * sqrt(vd)
  ub <- res.profile$b +  1.96 * sqrt(vd)
  
  ## Marginal AIC: OKAY to use if the main interest is to model fixed population effects with a reasonable correlation structur (Kneib & Greven (2010))
  mAIC <- 2 * res.profile$lk + 2 * (px + (length(mypar) - ns))  
  
  ## Conditional AIC: Vaida & Blanchard (2005) expression details in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2572765/pdf/nihms62635.pdf
  ## However a better approach should be in Kneib & Steven (2010) Appendix:B https://arxiv.org/pdf/1803.05664.pdf  
  ## assuming V and sigma^2 are (UN)KNOWN
  # if(pooled == T & common.s2 == T){  
  #   c2 <- diag(s2, ncol(V)) * V 
  #   
  #   #if(common.s2 == T){    
  #   #  c2 <- bdiag(lapply(seq_len(K), function(kk) diag(s2[1], ncol(V)) * V))  # 
  #   #} else if(common.s2 == F){
  #   #  c2 <- bdiag(lapply(seq_len(K), function(kk) diag(s2[kk], ncol(V)) * V)) # block diagonal 
  #   #}
  #   
  #   LLinv <- solve(c2, diag(1, ncol(c2)))
  #   c3 <- eigen(LLinv)
  #   Lambda <- c3$vectors %*% diag(sqrt(c3$values)) 
  #   
  #   c11 <- cbind(X, Z)
  #   c12 <- cbind(matrix(0, ncol = px, nrow = nrow(Lambda)), Lambda)
  #   # dim(c11); dim(c12)
  #   M <- rbind(c11, c12)
  #   
  #   MtM <- solve(t(M) %*% M)
  #   # H1 <- c11 %*% MtM %*% t(c11)
  #   # trace <- sum(diag(H1))
  #   trace <- sum(diag(MtM %*% t(c11) %*% c11))
  #   cAIC_vb <- 2 * res.profile$lk + 2 * trace 
  # } else{
  cAIC_vb = NULL
  # }
  
  ## Prediction 
  uihat <- as.matrix(do.call(rbind, lapply(seq_len(K), function(kk) { 
    
    ll <- length(which(id.site == id.site.uniq[kk]))
    dm <- matrix(1, nrow = ll, ncol = pz)
    
    sweep(dm, MARGIN = 2, res.profile$ui[[kk]], `*`)
    
  })))
  
  uihat_subj <- as.matrix(do.call(rbind, lapply(seq_len(K), function(kk) { 
    t(res.profile$ui[[kk]])
  })))
  
  if(pooled == T){
    Yhat <- X %*% as.matrix(res.profile$b) # population-level
    Yihat <- Yhat +  rowSums(Z * uihat) # subject-level
  } else {
    Yhat <- Yihat <- NULL 
  }

  return(list(b = res.profile$b, 
              b.sd = sqrt(vd),     # sd of fixed effect est
              wald = wald,   # Wald-test statistic
              lb = lb,       # lower-bound
              ub = ub,       # uppper-bound
              XvX = res.profile$allterms$bterm1, # X^{T}V^{-1}X
              ui = res.profile$ui, # BLUP of random effects
              uiM = uihat_subj, 
              varui = res.profile$varui,  # Variance (based on prediction error)
              varui_post = res.profile$varui_post, # posterior 
              Yhat = Yhat,         # population-level prediction (WOT random effects) XB
              Yihat = Yihat,       # subject-specific prediction XB + Zu
              mAIC = mAIC, 
              cAIC_vb = cAIC_vb, 
              V = V,
              s2 = s2,
              res = res, res.profile = res.profile))
}


#' @keywords internal
lmm.lrt <- function(Y = NULL, X = NULL, id.site = NULL, # Z = NULL, 
                    p.threshold = 0.05,
                    pooled = F, reml = F, 
                    common.s2 = T,      
                    SiXYZ = list(),
                    corstr = 'independence', # 'exchangeable', 'ar1', 'unstructured'),
                    # hessian = F,
                    verbose){
  ## model selection: select significant random effect of the covariates by LRT
  # random effect of covariate m at site i = u_im ~ N(0, vm), i=1...K, m=1...q
  # for each m, test H0: vm=0 vs H1: vm>0  
  # this is LRT on boundary, LR ~ 0.5*chisq(df=0) + 0.5*chisq(df=1) 
  
  if(pooled == T){
    N <- length(Y)
    id.site.uniq <- unique(id.site)
    px <- ncol(X)
    # pz <- ncol(Z)
    K <- length(id.site.uniq)
    SiXYZ <- list()
    for(ii in seq_along(id.site.uniq)){
      si <- id.site.uniq[ii]
      Xi <- X[id.site == si, ]
      # Zi <- Z[id.site == si, ]
      Yi <- Y[id.site == si]
      
      if(any(apply(Xi[,-1], 2, function(a)length(unique(a)))==1))
        warning(paste0('singular X in site #', ii, ' detected, REML maybe problematic!'))
      # if(any(apply(Zi[,-1], 2, function(a)length(unique(a)))==1))
      #  warning(paste0('singular X in site #', ii, ' detected, REML maybe problematic!'))
      
      SiXYZ[[si]]$SiX <- t(Xi) %*% Xi
      # SiXYZ[[si]]$SiXZ <- t(Xi) %*% Zi
      SiXYZ[[si]]$SiXY <- t(Xi) %*% Yi
      # SiXYZ[[si]]$SiZ  <- t(Zi) %*% Zi
      # SiXYZ[[si]]$SiZY <- t(Zi) %*% Yi
      SiXYZ[[si]]$SiY  <- sum(Yi ^ 2)
      SiXYZ[[si]]$ni <- sum(id.site == si)
    }
  }else{ 
    N <- sum(unlist(lapply(SiXYZ, function(a) a$ni)))
    id.site.uniq <- names(SiXYZ)
    px <- ncol(SiXYZ[[1]]$SiX)
    # pz <- ncol(SiXYZ[[1]]$SiXZ)
    K <- length(SiXYZ)
  }
  
  # LMM under H0 (fixed effect + random intercept), assume common residual var
  for(si in id.site.uniq){
    SiXYZ[[si]]$SiXZ <- as.matrix(SiXYZ[[si]]$SiX[,1])
    SiXYZ[[si]]$SiZ  <- as.matrix(SiXYZ[[si]]$SiX[1,1])
    SiXYZ[[si]]$SiZY <- as.matrix(SiXYZ[[si]]$SiXY[1,])
  }
  fit0 <- lmm.fit(SiXYZ = SiXYZ, pooled=pooled, common.s2 = common.s2, hessian=F, reml = reml)
  LR <- matrix(NA, px-1, 2)
  for(ii in 1:(px-1)){
    cat('covariate #', ii, '...\n')
    # LMM under H1 
    for(si in id.site.uniq){
      SiXYZ[[si]]$SiXZ <- as.matrix(SiXYZ[[si]]$SiX[,c(1,ii+1)])
      SiXYZ[[si]]$SiZ  <- as.matrix(SiXYZ[[si]]$SiX[c(1,ii+1),c(1,ii+1)])
      SiXYZ[[si]]$SiZY <- as.matrix(SiXYZ[[si]]$SiXY[c(1,ii+1),])
    }
    fit1i <- lmm.fit(SiXYZ = SiXYZ, pooled=pooled, common.s2 = common.s2, hessian=F, reml = reml)
    LR[ii,1] <- 2*(fit1i$res.profile$lp -fit0$res.profile$lp)
    LR[ii,2] <- 0.5 * (1-pchisq(LR[ii,1], df=1)) 
    if(LR[ii, 2]<p.threshold) cat('pval=', LR[ii, 2], ', R.E. needed \n')
  }    
  
  id.re <- which(LR[,2] < p.threshold)  # index for random effects
  for(si in id.site.uniq){
    SiXYZ[[si]]$SiXZ <- as.matrix(SiXYZ[[si]]$SiX[,c(1,id.re+1)])
    SiXYZ[[si]]$SiZ  <- as.matrix(SiXYZ[[si]]$SiX[c(1,id.re+1),c(1,id.re+1)])
    SiXYZ[[si]]$SiZY <- as.matrix(SiXYZ[[si]]$SiXY[c(1,id.re+1),])
  }
  fit.selected <- lmm.fit(SiXYZ = SiXYZ, pooled=pooled, common.s2 = common.s2, hessian=F, reml = reml)
  
  return(list(LR=LR, id.re=id.re, fit.selected=fit.selected))
}


#' @keywords internal
lmm.get.summary <- function(Y = NULL, X = NULL, Z = NULL, 
                            weights = NULL, id.site = NULL){
  ## get summary stats from each site for distributed lmm
  # 20201203: incorporate weight (wt) in LMM for distributed PQL, all the summary stats are Xi^TWiX_i, Xi^TWiZ_i, Xi^TWiY_i, etc
  
  if(is.null(weights)) weights <- rep(1, length(Y))
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  id.site <- as.character(id.site)
  id.site.uniq <- unique(id.site)
  px <- ncol(X)
  pz <- ncol(Z)
  
  SiXYZ <- list()
  for(ii in seq_along(id.site.uniq)){
    si = id.site.uniq[ii]
    wti = weights[id.site == si]
    Xi <- X[id.site == si, ]
    Zi <- Z[id.site == si, ]
    Yi <- Y[id.site == si]
    # if(any(apply(Xi[,-1], 2, function(a)length(unique(a)))==1)) 
    #   warning(paste0('singular X in site #', ii, ' detected!'))
    # if(any(apply(Zi[,-1], 2, function(a)length(unique(a)))==1)) 
    #   warning(paste0('singular Z in site #', ii, ' detected!'))
    
    SiX  = t(Xi*wti) %*% Xi
    SiXZ = t(Xi*wti) %*% Zi
    SiXY = t(Xi*wti) %*% Yi
    SiZ  = t(Zi*wti) %*% Zi
    SiZY = t(Zi*wti) %*% Yi
    SiY  = sum(Yi ^ 2 *wti)
    ni <- sum(id.site == si)
    SiXYZ[[si]] <- list(SiX  = SiX, SiXZ = SiXZ, SiXY = SiXY,
                        SiZ  = SiZ, SiZY = SiZY, SiY  = SiY, ni = ni)
  }
  
  return(SiXYZ)
}


#' @keywords internal
lmm.fit.partial <- function(id.re, SiXYZ, pooled=F, reml=T, hessian=T){
  for(si in names(SiXYZ)){
    SiXYZ[[si]]$SiZ = as.matrix(SiXYZ[[si]]$SiX[id.re, id.re])
    SiXYZ[[si]]$SiXZ = as.matrix(SiXYZ[[si]]$SiX[,id.re] ) 
    SiXYZ[[si]]$SiZY = as.matrix(SiXYZ[[si]]$SiXY[id.re])
  }
  fit1i <- lmm.fit(SiXYZ = SiXYZ, pooled=pooled, reml=reml, hessian=hessian)
  return(fit1i)
}


#' @keywords internal
lmm.get.summary.cut <- function(SiXYZ, site=NULL, x.fix=NULL, x.random=NULL){
  ## set DLMM summary stats by trimming site, x.fix and x.random
  
  if(is.null(x.fix)) site=1:length(SiXYZ)
  if(is.null(x.fix)) x.fix=1:nrow(SiXYZ[[1]]$SiX)
  if(is.null(x.random)) x.random=c(1) 
  
  if(is.character(site)) site = match(site, names(SiXYZ))
  if(is.character(x.fix)) x.fix = match(x.fix, colnames(SiXYZ[[1]]$SiX))
  if(is.character(x.random)) x.random = match(x.random, colnames(SiXYZ[[1]]$SiX))
  
  SiXYZ_ = SiXYZ[site]
  for(si in names(SiXYZ_)){
    SiXYZ_[[si]]$SiX =  SiXYZ[[si]]$SiX[x.fix, x.fix]
    SiXYZ_[[si]]$SiXY =  as.matrix(SiXYZ[[si]]$SiXY[x.fix])
    SiXYZ_[[si]]$SiZ = as.matrix(SiXYZ[[si]]$SiX[x.random, x.random])
    SiXYZ_[[si]]$SiXZ = as.matrix(SiXYZ[[si]]$SiX[,x.random] ) 
    SiXYZ_[[si]]$SiZY = as.matrix(SiXYZ[[si]]$SiXY[x.random])
  }
  return(SiXYZ_)
}


#' @keywords internal
dlmm <- function(SiXYZ, site=NULL, x.fix=NULL, x.random=NULL, 
                 select.re=NULL, select.re.pval = 0.05, # c('forward', 'univariate')
                 reml=T, hessian=T, table1=T, verbose=T){
  tab1 = NULL
  if(table1==TRUE){
    xnames <- colnames(SiXYZ[[1]]$SiX)
    px <- nrow(SiXYZ[[1]]$SiX)
    tab1.n <- matrix(unlist(lapply(SiXYZ, function(a) c(a$ni, a$SiX[-1,1], round(a$SiXY[1]/a$ni,1)))), nrow=px+1)
    tab1.p <- matrix(unlist(lapply(SiXYZ, function(a) round(c(a$ni/a$ni*100, a$SiX[-1,1]/a$ni*100, sqrt((a$SiY- a$SiXY[1]^2/a$ni)/(a$ni-1))),1))), nrow=px+1)
    tab1 = data.frame(matrix(paste0(tab1.n, ' (', tab1.p, ')'), nrow=px+1))
    names(tab1) <- names(SiXYZ)
    rownames(tab1) <- c('Total', xnames[-1], 'LOS')
  }
  
  if(is.null(site))  site <- c(1:length(SiXYZ))
  if(is.null(x.fix))  x.fix <- c(1:ncol(SiXYZ[[1]]$SiX))
  if(is.null(x.random))  x.random <- c(1)
  
  if(is.character(site))  site <- match(site, names(SiXYZ))
  if(is.character(x.fix))  x.fix <- match(x.fix, colnames(SiXYZ[[1]]$SiX))
  if(is.character(x.random))  x.random <- match(x.random, colnames(SiXYZ[[1]]$SiX))
  SiXYZ0 = lmm.get.summary.cut(SiXYZ, site=site, x.fix=x.fix, x.random=x.random)
  
  K <- length(SiXYZ0)
  xnames <- colnames(SiXYZ0[[1]]$SiX) 
  px <- length(xnames)
  nn <- unlist(lapply(SiXYZ0, function(a) a$ni))
  
  id.re <- x.random
  
  LR <- pval <- rep(NA, px-1)
  if(select.re=='univariate'){
    message('random-effects selection by univaraite LRT (based on random intercept)...')
    LR <- pval <- rep(NA, px-1)
    SiXYZ <- lmm.get.summary.cut(SiXYZ0, site=NULL, x.fix=NULL, x.random=c(1))
    fit0 = lmm.fit(SiXYZ = SiXYZ, pooled=F, reml=reml, hessian=hessian)
    for(ii in 1:(px-1)){  
      # cat('covariate #', ii, '=', xnames[ii+1], '...')
      # LMM under H1 
      SiXYZ = lmm.get.summary.cut(SiXYZ0, site=NULL, x.fix=NULL, x.random=c(1,ii+1))
      fit1i <- lmm.fit(SiXYZ = SiXYZ, pooled=F, reml=reml, hessian=hessian)
      # LR[ii,1] <- 2*(summary(fit1i)$logLik -  summary(fit0)$logLik)
      LR[ii] <- 2*(fit1i$res.profile$lk - fit0$res.profile$lk)
      pval[ii] <- 0.5 * (1-pchisq(LR[ii], df=1)) 
      # if(pval[ii]<select.re.pval) cat('pval=', LR[ii, 2], ', R.E. needed \n')
    } 
    LR <- data.frame(xnames=xnames[-1], x.add=2:px, LR=round(LR,1), pval=round(pval,3))
    # write.csv(LR, file='/Users/chl18019/Dropbox/R/DLMM/OHDSI/fig&tab/LRT_16_univ.csv')
    id.re <- c(1, which(pval<select.re.pval)+1)
  }else if(select.re=='forward'){
    message('random-effects selection by forward LRT (based on random intercept)...')
    id <- c(1) 
    pval1 = 0
    id1 = id
    LR1 = c(NA)
    while(pval1 < select.re.pval){
      SiXYZ = lmm.get.summary.cut(SiXYZ0, site=NULL, x.fix=NULL, x.random=id)
      fit0 <- lmm.fit(SiXYZ = SiXYZ, pooled=F, reml=reml, hessian=hessian)
      # fit0 = lmm.fit.partial(id, SiXYZ)
      id.add = setdiff(1:px, id)
      LR.add = rep(NA, px)
      for(idi in id.add){
        idt = sort(c(id, idi))
        SiXYZ = lmm.get.summary.cut(SiXYZ0, site=NULL, x.fix=NULL, x.random=idt)
        fit1i <- lmm.fit(SiXYZ = SiXYZ, pooled=F, reml=reml, hessian=hessian)
        # fit1i <- lmm.fit.partial(idt, SiXYZ)
        LR.add[idi] = 2*(fit1i$res.profile$lk - fit0$res.profile$lk)
      }
      pval1 = 0.5 * (1-pchisq(max(LR.add, na.rm=T), df=1))  
      id1 = c(id1, which.max(LR.add))
      LR1 = c(LR1, max(LR.add, na.rm=T))
      id = sort(id1)
      # cat(which.max(LR.add), ', ', xnames[which.max(LR.add)], ', ', max(LR.add, na.rm=T), ', ', pval, '\n')
    }
    pval <- 0.5 * (1-pchisq(LR1, df=1))  
    LR = data.frame(xnames=xnames[id1], x.add=id1, LR=round(LR1,1), pval=round(pval,3))
    id.re = c(1, sort(id1[pval < select.re.pval]) )
    # fit1.forward <- lmm.fit(SiXYZ = SiXYZ, pooled=F, reml=reml, hessian=hessian)
  } 
  
  SiXYZ <- lmm.get.summary.cut(SiXYZ0, site=NULL, x.fix=NULL, x.random=id.re)
  fit1 <- lmm.fit(SiXYZ = SiXYZ, pooled=F, reml=reml, hessian=hessian)
  
  return(list(SiXYZ=SiXYZ, LR=LR, table1=tab1, fit1=fit1))
  
}


#' @keywords internal
glmmDPQL.fit <- function(Y = NULL, X = NULL, Z = NULL, id.site = NULL, 
                         offset=NULL, weights = NULL, 
                         fixef.init = NULL, ranef.init=NULL,
                         pooled = F, reml = T, common.s2 = T,    # common residual var across sites
                         family, SiXYZ = list(),                 # IPD not available, use summ stats
                         corstr = 'independence',                # 'exchangeable', 'ar1', 'unstructured' 
                         hessian = F, 
                         niter = 10, tol=1e-06, verbose = T){
  ## distributed PQL (DPQL) via DLMM:
  ## MASS:glmmPQL essentially convert glmm estimation to iterative LMM
  ## DPQL utilize the lossless DLMM algorithm to replicate 
  ##      MASS:glmmPQL in a iterative and lossless way.
  ## notice the DLMM used is a weighted version, thus requires 
  ## t(Xi) %*% Wi %*% Xi, t(Xi) %*% Wi %*% Yi, t(Yi) %*% Wi %*% Yi
  ## from each site, where the weight Wi is updated in each iteration
  if (is.character(family)) 
    family <- get(family)
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if(is.null(weights)) weights=rep(1, length(Y))
  id.site.uniq <- unique(id.site)
  K <- length(id.site.uniq)
  qz <- ncol(Z)
  if(is.null(offset)) offset = rep(0, length(Y)) 
  
  # get init values via glm
  if(is.null(ranef.init)) ranef.init=matrix(0,K,qz)
  if(is.null(fixef.init)){  # use glm to get init working Y
    fit0 <- glm(Y ~ 0+X, offset=offset, family = family, weights=weights)
    w <- fit0$prior.weights
    eta <- fit0$linear.predictors
    zz <- eta + fit0$residuals  
    wz <- fit0$weights 
  } else{ 
    w <- weights
    eta = c(X%*%fixef.init )
    for(ii in seq_along(id.site.uniq))
      eta[id.site==id.site.uniq[ii]] <- eta[id.site==id.site.uniq[ii]] + as.matrix(Z[id.site==id.site.uniq[ii],]) %*% ranef.init[ii,]
    
    mu = family$linkinv(eta)
    mu.eta.val = family$mu.eta(eta)
    wz = w * mu.eta.val^2/family$variance(mu) 
    zz = eta + (Y-mu)/mu.eta.val - offset 
  }
  
  # formula.lmer = formula
  # formula.lmer[[2]] = quote(zz)
  u = ranef.init
  b = fixef.init
  bu =  c(b,u)
  
  AD.list <- list()
  
  for (i in seq_len(niter)) {
    if (verbose) message(gettextf("iteration %d", i), domain = NA) 
    # data$zz = zz
    # data$wz = wz
    # fit <- lmer(formula.lmer, data = data, REML = REML, weights = wz)
    
    AD.list[[i]] = lmm.get.summary(Y = zz, X = X, Z = Z, weights = wz, id.site = id.site)
    
    fit <- lmm.fit(Y = zz, X, Z, id.site, weights = wz, pooled = pooled, 
                   reml = reml, common.s2 = common.s2,  hessian = hessian, verbose=verbose)
    # etaold <- eta 
    eta <- X %*% fit$b  
    bu = cbind(bu, c(fit$b, unlist(fit$ui)))
    for(ii in seq_along(id.site.uniq))
      eta[id.site==id.site.uniq[ii]] <- eta[id.site==id.site.uniq[ii]] + as.matrix(Z[id.site==id.site.uniq[ii],]) %*% fit$ui[[ii]]
    
    # if (sum((eta - etaold)^2) < tol * sum(eta^2))  break
    dif = sum((bu[,i+1] - bu[,i])^2) / sum(bu[,i]^2) 
    if(verbose) message('diff=', dif)
    if (dif < tol)  break
    
    mu <- family$linkinv(eta)
    mu.eta.val <- family$mu.eta(eta)
    zz <- eta + (Y - mu)/mu.eta.val - offset
    wz <- w * mu.eta.val^2/family$variance(mu) 
  }
  fit$iter = i
  fit$bu = bu
  fit$AD.list = AD.list
  return(fit)
}
 







#' @useDynLib pda
#' @title A flexible version of MASS::glmmPQL
#' 
#' @usage myglmmPQL(formula.glm, formula, offset=NULL, family, data, fixef.init = NULL, 
#'                  weights=NULL, REML=T, niter=10, verbose=T)
#' @param formula.glm formula used to fit \code{glm} for initial fixed effects 
#' @param formula formula used to fit iterative \code{lmer} in PQL algorithm
#' @param offset \code{glm} offset term
#' @param family \code{glm} family
#' @param data \code{glm} data 
#' @param fixef.init initial fixed effects estimates, set to zeros if NULL
#' @param weights \code{glm} weights
#' @param REML \code{lmer} logical scalar - Should the estimates be chosen to optimize the REML criterion (as opposed to the log-likelihood)?
#' @param niter \code{glmmPQL} maximum number of iterations.
#' @param verbose \code{glmmPQL} logical: print out record of iterations?
#' @details Use lme4::lmer instead of nlme::varFixed in PQL iteration to allow REML
#'  
#' @return  An object wiht the same format as \code{lmer}.
#' @keywords internal
myglmmPQL <- function(formula.glm, formula, offset=NULL, family, data, 
                      fixef.init = NULL, weights=NULL, REML=T, niter = 10, verbose = T){
  if (is.character(family)) 
    family <- get(family)
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  if(is.null(offset)) offset <- rep(0, nrow(data))
  data$offset <- offset 
  if(is.null(weights)) weights <- rep(1, nrow(data))
  data$weights <- weights 
  
  
  # get init values via glm
  # if(is.null(ranef.init)) ranef.init=matrix(0,K,qz)
  if(is.null(fixef.init)){  # use glm to get init working Y
    fit0 <- glm(formula.glm, offset=offset, family = family, data=data, weights=weights)
    w <- fit0$prior.weights
    eta <- fit0$linear.predictors
    zz <- eta + fit0$residuals  
    wz <- fit0$weights
  } else{ 
    w <- weights
    X <- model.matrix(formula.glm, model.frame(formula.glm, data))
    Y <- as.numeric(model.response(model.frame(formula.glm, data)))
    eta <- c(X%*%fixef.init)
    # for(ii in seq_along(id.site.uniq))
    #   eta[id.site==id.site.uniq[ii]] <- eta[id.site==id.site.uniq[ii]] + as.matrix(Z[id.site==id.site.uniq[ii],]) %*% ranef.init[ii,]
    mu <- family$linkinv(eta)
    mu.eta.val <- family$mu.eta(eta)
    wz <- w * mu.eta.val^2/family$variance(mu) 
    zz <- eta + (Y-mu)/mu.eta.val - offset 
  }
  
  
  formula.lmer = formula
  formula.lmer[[2]] = quote(zz)
  for (i in seq_len(niter)) {
    if (verbose) message(gettextf("iteration %d", i), domain = NA)
    # fit <- eval(mcall)
    data$zz = zz
    data$wz = wz
    fit <- lme4::lmer(formula.lmer, data = data, REML = REML, weights = wz)
    etaold <- eta
    eta <- fitted(fit) + offset
    if (sum((eta - etaold)^2) < 1e-06 * sum(eta^2))  break
    mu <- family$linkinv(eta)
    mu.eta.val <- family$mu.eta(eta)
    zz <- eta + (Y - mu)/mu.eta.val - offset
    wz <- w * mu.eta.val^2/family$variance(mu) 
  }
  
  # attributes(fit$ logLik) <- NULL
  # fit$call <- Call
  # fit$family <- family
  # fit$logLik <- as.numeric(NA)
  # oldClass(fit) <- c("glmmPQL", oldClass(fit))
  return(fit)
}





# fit.pool <- myglmmPQL(death~age+sex+lab, death~age+sex+lab+(1|site), data=covid, 
#                       fixef.init=rep(0,4), family='binomial')
# fit.pool <- glmmDPQL.fit(Y=covid$death, X=cbind(1,as.matrix(covid[,2:4])), Z=matrix(1,nrow(covid),1), id.site = covid$site,
#                          fixef.init=rep(0,4), family='binomial', pooled = T, niter = control$maxround)
# cbind(bu.pool=c(fit.pool$b, unlist(fit.pool$ui)),
#       bu.dpql=c(fit.dpql$bhat,fit.dpql$uhat),
#       sd.pool=c(fit.pool$b.sd, sqrt(unlist(fit.pool$varui_post))),
#       sd.dpql=c(fit.dpql$sebhat, fit.dpql$seuhat))
# c(fit.pool$V, fit.dpql$Vhat)