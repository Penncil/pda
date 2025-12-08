
## One-shot Distributed Algorithms for Cox regression with Time-varying coefficients (ODACT)
# C. Jason Liang, Chongliang Luo

# Local constant partial likelihood
# Cai and Sun (2003) propose maximizing a local constant partial likelihood to estimate the time-dependent log-hazard ratio β(t) from a Cox regression model with time-dependent β(t):
# 
#   λ(t|x)=λ0(t)exp(β(t)x)
# Note that the “local” in “local constant partial likelihood” refers to “local in time” rather than a “local site” in the context of a distributed algorithm.
# 
# The local constant log partial likelihood is defined as
# L(β,t)=(nhn)−1∑i=1nK(yi−thn)δi[βTxi−log{∑j∈R(ti)exp(βTxj)}],
# where solving for β provides an estimate of β(t). The K(⋅) function is a kernel, such as the Epanechnikov kernel, centered around t with data-adaptive bandwidth hn.
# 
# Note that L(β,t) is similar to the usual partial log likelihood L(β), with the main difference being that each likelihood contribution is now weighted by a kernel centered around t. It then follows that the gradient and hessian of L(β,t) can be written as locally weighted versions of ∇L(β) and ∇2L(β):
#   ∇L(β,t)∇2L(β,t)=(nhn)−1∑i=1nK(yi−thn)δi{xi−∑j∈R(ti)exp(βTxj)xj∑j∈R(ti)exp(βTxj)}=(nhn)−1∑i=1nK(yi−thn)δi⎡⎣⎢⎢{∑j∈R(ti)exp(βTxj)xj}⨂2−∑j∈R(ti)exp(βTxj)∑j∈R(ti)exp(βTxj)x⨂2j{∑j∈R(ti)exp(βTxj)}2⎤⎦⎥⎥
# ## Functions for local in time estimation:


require(RcppArmadillo)
require(Rcpp)
# sourceCpp('pda/src/odact.cpp')

## DO NOT use order.time=T, always order time before calling 
## llpl llplg llplh sllpl sllplg D12

# log-lik of one site
#' @keywords internal
llpl <- function(beta, times, status, covars, tt=0.5, h=500,
                 order.time=F) {  
  # locally constant log partial likelihood
  if(order.time==T){
    times = times[order(times)]
    status = status[order(times)]
    covars = covars[order(times), ,drop=FALSE]
  }
  
  nn <- 1 - ((times - tt)/h)^2 >= 0
  K <- (1 - ((times - tt)/h)^2) * 0.75/h
  N <- length(times)
  if (is.null(ncol(covars))) 
    lp <- covars * beta
  else lp <- covars %*% beta
  
  # -sum((K * status * (lp - log(cumsum(exp(lp)[N:1])[N:1])))[nn])
  ## cum_sum is rcpp function, ~ 20 times faster
  -sum((K * status * (lp - log(cum_sum(exp(lp), reversely = T))))[nn])
}

# gradient of ll of one site
#' @keywords internal
llplg <- function(beta, times, status, covars, tt=0.5, h=500,
                  order.time=F) {
  # locally constant log partial likelihood gradient
  if(order.time==T){
    times = times[order(times)]
    status = status[order(times)]
    covars = covars[order(times), ,drop=FALSE]
  }
  
  nn <- 1 - ((times - tt)/h)^2 >= 0
  K <- (1 - ((times - tt)/h)^2) * 0.75/h
  N <- length(times)
  if (is.null(ncol(covars))) 
    lp <- covars * beta
  else lp <- covars %*% beta
  elp <- c(exp(lp))
  # num <- apply(elp*covars, 2, function(x) rev(cumsum(rev(x))))
  # den <- cumsum(elp[N:1])[N:1]   # %o% rep(1, ncol(covars))
  ## cum_sum_col is rcpp function, ~ 20 times faster
  num <- cum_sum_cols(elp*covars, reversely = T)
  den <- c(cum_sum(elp, reversely = T))
  
  -colSums(((covars - num/den) * (K * status))[nn, , drop = FALSE])
}

# hessian of ll of one site
#' @keywords internal
llplh <- function(beta, times, status, covars, tt=0.5, h=500,
                  order.time=F) {
  # locally constant log partial likelihood gradient
  if(order.time==T){
    times = times[order(times)]
    status = status[order(times)]
    covars = covars[order(times), ,drop=FALSE]
  }
  
  nn <- 1 - ((times - tt)/h)^2 >= 0
  K <- (1 - ((times - tt)/h)^2) * 0.75/h
  N <- length(times)
  px = ncol(covars)
  if (is.null(ncol(covars))){
    lp <- covars * beta
  }else{
    lp <- covars %*% beta
  }
  
  elp <- c(exp(lp))
  
  # den <- ((cumsum(elp[N:1])[N:1])^2)
  # num1 <- apply(elp*covars, 2, function(x) rev(cumsum(rev(x))))
  # num1 = t(apply(num1, 1, function(x) c(x %o% x ) ))   # matrix(, nrow=N)
  # num2 = elp * t(apply(covars, 1, function(x) c(x %o% x) ) )
  # num2 <- sqrt(den) * apply(num2, 2, function(x) rev(cumsum(rev(x))) )
  
  ## cum_sum_col Xotimes2 are rcpp function, ~ 10 times faster
  den <- c(cum_sum(elp, reversely = T)^2)
  num1 <- cum_sum_cols(elp*covars, reversely = T)
  num1 <- Xotimes2(num1)
  num2 <- elp * Xotimes2(covars)
  num2 <- sqrt(den) * cum_sum_cols(num2, reversely = T)
  
  hh = -colSums((((num1-num2)/den) * (K * status))[nn, , drop = FALSE])
  hh = matrix(hh, px, px)
  
  return(hh)
}

# fit Cox with beta(t) for one site
#' @keywords internal
seq_fit_list <- function(da, fn=llpl, times=seq(0,1,0.1), h=0.1, betabar, ...){
  # wrapper function to make it easier to call optim() multiple times
  # da: data table, columns: time, status, X's
  # return: list of optim output
  lapply(as.list(times), function(x)
    optim(par=betabar, fn=fn, method="Nelder-Mead", times=da[,1], status=da[,2], 
          covars=as.matrix(da[,-c(1:2),drop=FALSE]), tt=x, h=h, ...) )
}


# log-lik of site-stratified Cox with beta(t)
#' @keywords internal
llpl_st <- function(beta, times, status, covars, site, tt=0.5, h=500,
                    order.time=F){
  site.uni = unique(site)
  sum(sapply(site.uni, function(ss) llpl(beta, times[site==ss], status[site==ss], covars[site==ss,],
                                         tt=tt, h=h, order.time=order.time) ) )
}

# fit site-stratified Cox with beta(t), list of optim outputs  
#' @keywords internal
seq_fit_st_list <- function(da, site, fn=llpl_st, times=seq(0,1,0.1), h=0.1, betabar, ...){
  # wrapper function to make it easier to call optim() multiple times
  da = data.frame(da)
  lapply(as.list(times), function(x)
    optim(par=betabar, fn=fn, method="Nelder-Mead", times=da[,1], status=da[,2], covars=as.matrix(da[,-c(1:2),drop=FALSE]),
          site=site, tt=x, h=h, ...) )
}

# fit site-stratified Cox with beta(t), as matrix format  
#' @keywords internal
mycoxph_bt <- function(da, site, fn=llpl_st, times=seq(0,1,0.1), h=0.1, betabar, ...){
  px = ncol(da) - 2
  fit.pool = seq_fit_st_list(da, site, llpl_st, times, h, betabar, ...)
  b.pool = matrix(unlist(lapply(fit.pool, function(a) a$par) ), nrow=px) 
  se.pool = sqrt(matrix(unlist(sapply(fit.pool, function(a) diag(solve(a$hessian)*0.6/h) )), nrow=px) )
  return(list(b.pool=round(b.pool,4), se.pool= round(se.pool,4)))
}
  
# # surrogate ll 
# sllpl <- function(beta, times, status, covars, betabar, ind1, order=1, tt=0.5, h=.25,
#                   order.time=F) {
#   # surrogate locally constant log partial likelihood
#   if(order.time==T){
#     times = times[order(times)]
#     status = status[order(times)]
#     covars = covars[order(times), ,drop=FALSE]
#   }
#   
#   if(order==1){
#     # first order version (set 'order=1')
#     llpl(beta, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1]) +
#       sum(beta * (llplg(betabar, times=times, status=status, covars=covars, tt=tt, h=h)/length(times) -
#                     llplg(betabar, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1])))
#   } else{
#     # second order version (set 'order=2')
#     llpl(beta, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1]) +
#       sum(beta * (llplg(betabar, times=times, status=status, covars=covars, tt=tt, h=h)/length(times) - 
#                     llplg(betabar, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1]))) +
#       0.5 * t(beta - betabar) %*% (llplh(betabar, times=times, status=status, covars=covars, tt=tt, h=h)/length(times) - 
#                                      llplh(betabar, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1])) %*% (beta - betabar)
#     
#   }
# }
# 
# # gradient of sll
# sllplg <- function(beta, times, status, covars, betabar, ind1, order=1, tt=0.5, h=.25,
#                    order.time=F) {
#   # surrogate locally constant log partial likelihood
#   if(order.time==T){
#     times = times[order(times)]
#     status = status[order(times)]
#     covars = covars[order(times), ,drop=FALSE]
#   }
#   
#   if(order==1){
#     # first order version (set 'order=1')
#     llplg(beta, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1]) +
#       (llplg(betabar, times=times, status=status, covars=covars, tt=tt, h=h)/length(times) -
#          llplg(betabar, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1])) 
#   } else{
#     # second order version (set 'order=2')
#     llplg(beta, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1]) +
#       (llplg(betabar, times=times, status=status, covars=covars, tt=tt, h=h)/length(times) - 
#          llplg(betabar, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1])) +
#       0.5 * (llplh(betabar, times=times, status=status, covars=covars, tt=tt, h=h)/length(times) - 
#                llplh(betabar, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1])) %*% (beta - betabar)
#     
#   }
# }


# 
# # gradient of llpl_st
# llplg_st <- function(beta, times, status, covars, site, tt=0.5, h=500,
#                      order.time=F) {
#   site.uni = unique(site)
#   rowSums(sapply(site.uni, function(ss) llplg(beta, times[site==ss], status[site==ss], covars[site==ss,], tt=tt, h=h,
#                                               order.time=order.time) ) )
# }
# # llplg_st(cox_st$coef, dmat[,1], dmat[,2], dmat[,-c(1,2)], dm0$site, 1,2)
# 
# # hessian of ll
# llplh_st <- function(beta, times, status, covars, site, tt=0.5, h=500,
#                      order.time=F) {
#   site.uni = unique(site)
#   px = length(beta)
#   tmp = sapply(site.uni, function(ss) llplh(beta, times[site==ss], status[site==ss], covars[site==ss,], tt=tt, h=h,
#                                             order.time=order.time) )  
#   matrix(rowSums(tmp), px, px)
# } 
# # llplh_st(cox_st$coef[1:2], dmat[,1], dmat[,2], dmat[,3:4], dm0$site, 1,2)
# 
# # surrogate ll 
# sllpl_st <- function(beta, times, status, covars, site, betabar, ind1, order=1, tt=0.5, h=.25,
#                      order.time=F) {
#   # surrogate locally constant log partial likelihood
#   if(order.time==T){
#     times = times[order(times)]
#     status = status[order(times)]
#     covars = covars[order(times), ,drop=FALSE]
#   }
#   
#   if(order==1){
#     # first order version (set 'order=1')
#     llpl(beta, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1]) +
#       sum(beta * (llplg_st(betabar, times=times, status=status, covars=covars, site=site, tt=tt, h=h)/length(times) -
#                     llplg(betabar, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1])))
#   } else{
#     # second order version (set 'order=2')
#     llpl(beta, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1]) +
#       sum(beta * (llplg_st(betabar, times=times, status=status, covars=covars, site=site, tt=tt, h=h)/length(times) - 
#                     llplg(betabar, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1]))) +
#       0.5 * t(beta - betabar) %*% (llplh_st(betabar, times=times, status=status, covars=covars, site=site, tt=tt, h=h)/length(times) - 
#                                      llplh(betabar, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1])) %*% (beta - betabar)
#     
#   }
# }
# 
# # gradient of sll
# sllplg_st <- function(beta, times, status, covars, site, betabar, ind1, order=1, tt=0.5, h=.25,
#                       order.time=F) {
#   # surrogate locally constant log partial likelihood
#   if(order.time==T){
#     times = times[order(times)]
#     status = status[order(times)]
#     covars = covars[order(times), ,drop=FALSE]
#   }
#   
#   if(order==1){
#     # first order version (set 'order=1')
#     llplg(beta, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1]) +
#       (llplg_st(betabar, times=times, status=status, covars=covars, site=site, tt=tt, h=h)/length(times) -
#          llplg(betabar, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1])) 
#   } else{
#     # second order version (set 'order=2')
#     llplg(beta, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1]) +
#       (llplg_st(betabar, times=times, status=status, covars=covars, site=site, tt=tt, h=h)/length(times) - 
#          llplg(betabar, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1])) +
#       0.5 * (llplh_st(betabar, times=times, status=status, covars=covars, site=site, tt=tt, h=h)/length(times) - 
#                llplh(betabar, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1])) %*% (beta - betabar)
#     
#   }
# }
 

 

# # 1st and 2nd order derivatives of llpl, using local or pooled data
# D12 <- function(times, status, covars, betabar, ind1,  tt=0.5, h=.25,
#                 order.time=F) {
#   # output derivatives at betabar
#   if(order.time==T){
#     times = times[order(times)]
#     status = status[order(times)]
#     covars = covars[order(times), ,drop=FALSE]
#   }
#   
#   # second order version (set 'order=2')
#   D1.N = llplg(betabar, times=times, status=status, covars=covars, tt=tt, h=h)/length(times)  
#   D1.1 = llplg(betabar, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1]) 
#   D2.N = llplh(betabar, times=times, status=status, covars=covars, tt=tt, h=h)/length(times)  
#   D2.1 = llplh(betabar, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1]) 
#   
#   return(list(D1.N=D1.N, D1.1=D1.1, D2.N=D2.N, D2.1=D2.1))
# }

# seq_fit <- function(da, fn=llpl, times=seq(0,1,0.1), h=0.1, betabar, ...){
#   # wrapper function to make it easier to call optim() multiple times
#   # da: data table, columns: time, status, X's
#   # return: beta matrix of beta * times
#   apply(as.matrix(times), 1, function(x)
#     optim(par=betabar, fn=fn, method="Nelder-Mead", times=da[,1], status=da[,2], 
#           covars=da[,-c(1:2),drop=FALSE], tt=x, h=h, ...)$par)
# }







# # generate time-varying cox data
# # lambda(t|m) = exp(tm+m) w/ baseline haz = 1
# # M ~ Unif(0,1)
# Stm <- function(t,m,c0=-1){
#   #S(t|m)
#   #exp((1/m) - (exp(t*m)/m))
#   # exp(exp(m)/m - exp(t*m + m)/m)
#   exp(exp(m*c0)/m*(1- exp(t*m)) )
# } 
# 
# # generate random times given m
# Stm_inv <- function(zo, m,c0=-1){
#   # helper function for generating random times given m
#   #zo for random number from Unif(0, 1)
#   #  log(1 - m * log(zo))/m
#   (log(exp(m*c0) - m*log(zo)) - m*c0)/m
# } 
# 
# ftm <- function(t,m){
#   # spot check to make sure random numbers are valid
#   # f(t, m)
#   #  if (m < 1 | m > 3) return(0)
#   #  (1/2) * exp((1/m) - (exp(t*m)/m) + t*m)
#   #(1/2) * exp(exp(m)/m - exp(t*m + m)/m + t*m + m)
#   exp(exp(m)/m - exp(t*m + m)/m + t*m + m)
# }
# 
# ft <- function(t){
#   # f(t), numerically integrated
#   # again, just spot checking
#   ms <- seq(0,1,.001)
#   sum(ftm(rep(t, length(ms)), ms))*.001
# }
 