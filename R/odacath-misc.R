#Heterogeneous ODAM Code
#Expit
#Needs to use R 3.5.1 or higher
#UPDATE: 10/12/2022 We now use the faster version of MLR log-likelihood to speed up calculation times in code 
#UPDATE2: 10/18/2022 Instead of using newton raphson to find roots of score equation with multiroot, we instead aim to minimize the L2 norm of S(theta) and estimate
#coefficients with optim
#UPDATE3: 12/16/2022 Included a new meta analysis meta_analysis_v3 that calculates the site-specific intercepts on the original scale
# library(numDeriv)
# library(rootSolve)
# library(minqa)
# library(ordinal)
# library(plyr)
## library(brglm2)
## library(UPG)
#' @import minqa ordinal plyr
#' @keywords internal
expit=function(a) 1/(1+exp(-a))

#Likelihood-Used for estimation of parameters within each site (All parameters are estimated
#x and y (within site)
#' @keywords internal
Lik2=function(g,x,y,model=NULL){
  #Model Types=mlr (Multinomial Logistic Regression), polr (Proportional Odds Logistic Regression), lr (Logistic Regression)
  x <- as.matrix(x)
  n <- nrow(x)
  
  px <- ncol(x)
  ind_px=1:px
  q <- length(unique(y)) - 1
  ind_q <- seq_len(q) 
  
  if(model=="mlr"){
    x = cbind(1,x)
    if(q <= 1L) stop("response must have 3 or more levels")
    bm = cbind(matrix(g, px+1, q),0)
    Xb = x %*% bm
    y_cats=sort(unique(y))
    li=lapply(y_cats,function(l){ #Break up calculation by outcome category
      yindices=which(y==l)
      Xby=x[yindices,] %*% bm[,l]
      Xb_sub=Xb[yindices,]
      if(length(yindices)==1){
        return(Xby-log(sum(exp(Xb_sub))))  
      }else{
        return(Xby-log(rowSums(exp(Xb_sub))))
      }})
    logL =  - sum(unlist(li))
  }else if(model=="polr"){
    if(q <= 1L) stop("response must have 3 or more levels")
    theta=g[px + ind_q]
    gamm=c(-Inf, cumsum(c(theta[1L], exp(theta[-1L]))), Inf)
    # eta <- offset
    # if(px ) eta <- eta + drop(x %*% g[ind_px])
    if(px ) ed=drop(x %*% g[ind_px])
    z1=pmin(100, gamm[y+1L] - ed)
    z2=pmax(-100, gamm[y] - ed)
    pr=expit(z1) - expit(z2)
    # logL <- if (all(pr > 0)) -sum(wt * log(pr)) else Inf
    logL=if (all(pr > 0)) -sum(log(pr)) else Inf
  }else{
    x = cbind(1,x)
    logL=sum(-1*(y * (x%*%t(t(g))) - log(1 + exp(x%*%t(t(g))))))
  }
  return(logL/n)
}

#Vector of parameter assignments given beta and eta,for use in all model types. g is a general vector of parameters to be able
#to run code for all model types
#' @keywords internal
g_assignment=function(beta,eta,eta_indices){
  if(is.null(eta)==FALSE){ #Is there site specific parameters (Heterogeneity?)
    l_beta=length(beta)
    l_eta=length(eta)
    l_tot=l_beta+l_eta
    g=rep(0,times=l_tot)
    beta_indices=seq(1:l_tot)
    beta_indices=beta_indices[-eta_indices]
    g[beta_indices]=beta
    g[eta_indices]=eta
  }else{
    g=beta
  }
  return(list(g=g,beta_indices=beta_indices,eta_indices=eta_indices))
}

#Likelihood-Used for estimation of efficient score function (Only beta is estimated) (Also can be used to update eta given a beta when T>1)
#x and y (within site)
#' @keywords internal
Lik=function(beta,eta=NULL,x,y,model=NULL,eta_indices=NULL){
  #Model Types=mlr (Multinomial Logistic Regression), polr (Proportional Odds Logistic Regression), lr (Logistic Regression)
  #eta_indices=indices of g which are considered eta parameters (site specific)
  x <- as.matrix(x)
  n <- nrow(x)
  
  px <- ncol(x)
  ind_px=1:px
  q <- length(unique(y)) - 1
  ind_q <- seq_len(q) 
  
  
  #Assignment of beta and eta to parameter vector
  g_list=g_assignment(beta=beta,eta=eta,eta_indices=eta_indices)
  g=g_list$g
  
  if(model=="mlr"){
    x = cbind(1,x)
    if(q <= 1L) stop("response must have 3 or more levels")
    bm = cbind(matrix(g, px+1, q),0)
    Xb = x %*% bm
    y_cats=sort(unique(y))
    li=lapply(y_cats,function(l){ #Break up calculation by outcome category
      yindices=which(y==l)
      Xby=x[yindices,] %*% bm[,l]
      Xb_sub=Xb[yindices,]
      if(length(yindices)==1){
        return(Xby-log(sum(exp(Xb_sub))))  
      }else{
        return(Xby-log(rowSums(exp(Xb_sub))))
      }})
    logL =  - sum(unlist(li))
  }else if(model=="polr"){
    if(q <= 1L) stop("response must have 3 or more levels")
    theta=g[px + ind_q]
    gamm=c(-Inf, cumsum(c(theta[1L], exp(theta[-1L]))), Inf)
    # eta <- offset
    # if(px ) eta <- eta + drop(x %*% g[ind_px])
    if(px ) ed=drop(x %*% g[ind_px])
    z1=pmin(100, gamm[y+1L] - ed)
    z2=pmax(-100, gamm[y] - ed)
    pr=expit(z1) - expit(z2)
    # logL <- if (all(pr > 0)) -sum(wt * log(pr)) else Inf
    logL=if (all(pr > 0)) -sum(log(pr)) else Inf
  }else{
    x = cbind(1,x)
    logL=sum(-1*(y * (x%*%t(t(g))) - log(1 + exp(x%*%t(t(g))))))
  }
  return(logL/n)
}

#pdfs
#x and y (within site)
#' @keywords internal
pdf=function(beta,eta=NULL,x,y,model=NULL,eta_indices=NULL){
  x <- as.matrix(x)
  n <- nrow(x)
  px <- ncol(x)
  ind_px=1:px
  q <- length(unique(y)) - 1
  ind_q <- seq_len(q) 
  
  #Assignment of beta and eta to parameter vector
  g_list=g_assignment(beta=beta,eta=eta,eta_indices=eta_indices)
  g=g_list$g
  
  if(model=='mlr'){
    x = cbind(1,x)
    bm = cbind(matrix(g, px+1, q),0)
    Xb = x %*% bm
    y_cats=sort(unique(y))
    pdf=vector(length=nrow(x))
    for(l in y_cats)
    {
      yindices=which(y==l) 
      Xb_sub=Xb[yindices,]
      Xby=Xb_sub[,l]
      if(length(yindices)==1){
        pdf[yindices]=exp(Xby)/sum(exp(Xb_sub))   
      }else{
        pdf[yindices]=exp(Xby)/rowSums(exp(Xb_sub))
      }}
  }else if (model=='polr'){
    theta=g[px + ind_q]
    gamm=c(-Inf, cumsum(c(theta[1L], exp(theta[-1L]))), Inf)
    # eta <- offset
    # if(px ) eta <- eta + drop(x %*% g[ind_px])
    if(px ) ed=drop(x %*% g[ind_px])
    z1=pmin(100, gamm[y+1L] - ed)
    z2=pmax(-100, gamm[y] - ed)
    pdf=expit(z1) - expit(z2)
  }else{
    x = cbind(1,x)
    pdf=as.vector((expit(x%*%t(t(g))))^y*(1-expit(x%*%t(t(g))))^(1-y))
  }
  return(pdf)
}

#Likelihood-Tilting Ratio
#x and y (within site)
#' @keywords internal
Lik_tilt=function(g,x,y,model=NULL,tilt_ratio){
  #Model Types=mlr (Multinomial Logistic Regression), polr (Proportional Odds Logistic Regression), lr (Logistic Regression)
  #eta_indices=indices of g which are considered eta parameters (site specific)
  x <- as.matrix(x)
  n <- nrow(x)
  px <- ncol(x)
  ind_px=1:px
  q <- length(unique(y)) - 1
  ind_q <- seq_len(q) 
  
  
  if(model=="mlr"){
    x = cbind(1,x)
    if(q <= 1L) stop("response must have 3 or more levels")
    bm = cbind(matrix(g, px+1, q),0)
    Xb = x %*% bm
    y_cats=sort(unique(y))
    li=lapply(y_cats,function(l){ #Break up calculation by outcome category
      yindices=which(y==l)
      Xby=x[yindices,] %*% bm[,l]
      tilt_ratio_y=tilt_ratio[yindices]
      Xb_sub=Xb[yindices,]
      if(length(yindices)==1){
        return((Xby-log(sum(exp(Xb_sub))))*tilt_ratio_y)  
      }else{
        return((Xby-log(rowSums(exp(Xb_sub))))*tilt_ratio_y)
      }})
    logL =  - sum(unlist(li))
  }else if(model=="polr"){
    if(q <= 1L) stop("response must have 3 or more levels")
    theta=g[px + ind_q]
    gamm=c(-Inf, cumsum(c(theta[1L], exp(theta[-1L]))), Inf)
    # eta <- offset
    # if(px ) eta <- eta + drop(x %*% g[ind_px])
    if(px ) ed=drop(x %*% g[ind_px])
    z1=pmin(100, gamm[y+1L] - ed)
    z2=pmax(-100, gamm[y] - ed)
    pr=expit(z1) - expit(z2)
    # logL <- if (all(pr > 0)) -sum(wt * log(pr)) else Inf
    li=log(pr)*tilt_ratio
    logL=if (all(pr > 0)) -sum(li) else Inf
  }else{
    x = cbind(1,x)
    li=(y * (x%*%t(t(g))) - log(1 + exp(x%*%t(t(g)))))*tilt_ratio
    logL=-sum(li)
  }
  return(logL/n)
}

#Likelihood-Tilting Ratio
#x and y (within site)
#' @keywords internal
Lik_tilt_es=function(g,x,y,model=NULL,tilt_ratio){ 
  #Model Types=mlr (Multinomial Logistic Regression), polr (Proportional Odds Logistic Regression), lr (Logistic Regression)
  #eta_indices=indices of g which are considered eta parameters (site specific)
  x <- as.matrix(x)
  n <- nrow(x)
  
  px <- ncol(x)
  ind_px=1:px
  q <- length(unique(y)) - 1
  ind_q <- seq_len(q) 
  
  
  if(model=="mlr"){
    x = cbind(1,x)
    if(q <= 1L) stop("response must have 3 or more levels")
    bm = cbind(matrix(g, px+1, q),0)
    Xb = x %*% bm
    y_cats=sort(unique(y))
    li=lapply(y_cats,function(l){ #Break up calculation by outcome category
      yindices=which(y==l)
      Xby=x[yindices,] %*% bm[,l]
      tilt_ratio_y=tilt_ratio[yindices]
      Xb_sub=Xb[yindices,]
      if(length(yindices)==1){
        return((Xby-log(sum(exp(Xb_sub))))*tilt_ratio_y)  
      }else{
        return((Xby-log(rowSums(exp(Xb_sub))))*tilt_ratio_y)
      }})
    logL =  - sum(unlist(li))
  }else if(model=="polr"){
    if(q <= 1L) stop("response must have 3 or more levels")
    theta=g[px + ind_q]
    gamm=c(-Inf, cumsum(c(theta[1L], exp(theta[-1L]))), Inf)
    # eta <- offset
    # if(px ) eta <- eta + drop(x %*% g[ind_px])
    if(px ) ed=drop(x %*% g[ind_px])
    z1=pmin(100, gamm[y+1L] - ed)
    z2=pmax(-100, gamm[y] - ed)
    pr=expit(z1) - expit(z2)
    # logL <- if (all(pr > 0)) -sum(wt * log(pr)) else Inf
    li=log(pr)*tilt_ratio
    logL=if (all(pr > 0)) -sum(li) else Inf
  }else{
    x = cbind(1,x)
    li=(y * (x%*%t(t(g))) - log(1 + exp(x%*%t(t(g)))))*tilt_ratio
    logL=-sum(li)
  }
  
  return(logL)
}


#x and y (within site)
#' @keywords internal
tilt_I=function(betabar,etaj,eta1,x,y,eta_indices,model){ #Calculating the tilted information matrix for all individuals contributions at site j
  
  pdf_j=pdf(beta=betabar,eta=etaj,x=x,y=y,model=model,eta_indices=eta_indices)
  pdf_1=pdf(beta=betabar,eta=eta1,x=x,y=y,model=model,eta_indices=eta_indices)
  tilt_ratio = pdf_j/pdf_1
  
  #Assignment of beta and eta to parameter vector
  g_list=g_assignment(beta=betabar,eta=etaj,eta_indices=eta_indices)
  g=g_list$g
  
  ll=function(g){Lik_tilt(g=g,x=x,y=y,model=model,tilt_ratio)}
  tilde_I = NULL
  tilde_I=numDeriv::hessian(ll, g)
  return(tilde_I)
}

#x and y (within site)
#' @keywords internal
regular_S=function(beta,eta, x, y,eta_indices,model){ #Score for logistic regression for current model parameterization, used to calculate fisher information within each site (No tilting ratio)
  #Assignment of beta and eta to parameter vector
  g_list=g_assignment(beta=beta,eta=eta,eta_indices=eta_indices)
  g=g_list$g
  n=nrow(x)
  ll=function(g){Lik2(g=g,x=x,y=y,model=model)}
  S = NULL
  S=numDeriv::grad(ll, g)
  return(matrix(-1*n*S,ncol=1,nrow=length(g))) #remember log likelihood is divided by n and multiplied by -1
}
#x and y (within site)
#' @keywords internal
regular_I=function(beta,eta, x, y,eta_indices,model){ #Fisher information for logistic regression for current model parameterization, used to calculate fisher information within each site (No tilting ratio)
  #Assignment of beta and eta to parameter vector
  g_list=g_assignment(beta=beta,eta=eta,eta_indices=eta_indices)
  g=g_list$g
  
  n=nrow(x)
  ll=function(g){Lik2(g=g,x=x,y=y,model=model)}
  I = NULL
  I=numDeriv::hessian(ll, g)
  return(n*I) #remember log likelihood is divided by n and multiplied by -1: Information is negative inverse hessian
}
#x and y (within site)
#' @keywords internal
regular_es = function(beta,eta,x,y,eta_indices,model){ #Efficient Score Function that is calculated at each site j, S_j(\bar{\beta},\bar{\gamma}_j), (Calculates Fisher Information and partitions out into blocks)
  l_beta = length(beta)
  l_eta=length(eta)
  l_tot=l_beta+l_eta
  n=nrow(x)
  g_list=g_assignment(beta=beta,eta=eta,eta_indices=eta_indices)
  beta_indices=g_list$beta_indices
  
  #Score Calculation
  S=regular_S(beta=beta,eta=eta, x=x, y=y,eta_indices=eta_indices,model=model)/n
  
  #Hessian Calculation
  I=regular_I(beta=beta,eta=eta, x=x, y=y,eta_indices=eta_indices,model=model)/n # (dividing by n, the likelihood is already negative)
  
  I_ee = I[eta_indices,eta_indices]
  I_be = I[beta_indices,eta_indices]
  
  #Score 
  S_b = S[beta_indices,] #First order gradient for shared beta
  S_e = S[eta_indices,] #First order gradient for gamma and eta
  return((S_b-I_be%*%solve(I_ee)%*%S_e))
}

#x and y within site
#' @keywords internal
tilt_es = function(beta,betabar,etaj,eta1,x,y,eta_indices,model){ #Calculates the tilted score function at the local site for site j #Gamma1 (Site j), Gamma2 (Site 1)
  pdf_j=pdf(beta=betabar,eta=etaj,x=x,y=y,model=model,eta_indices=eta_indices)
  pdf_1=pdf(beta=betabar,eta=eta1,x=x,y=y,model=model,eta_indices=eta_indices)
  tilt_ratio = pdf_j/pdf_1
  n=nrow(x)
  
  #Assignment of beta and etaj to parameter vector
  ##Betabar and etaj
  g_betabar_etaj_list=g_assignment(beta=betabar,eta=etaj,eta_indices=eta_indices)
  g_betabar_etaj=g_betabar_etaj_list$g
  beta_indices=g_betabar_etaj_list$beta_indices
  ##Beta and eta
  g_beta_etaj_list=g_assignment(beta=beta,eta=etaj,eta_indices=eta_indices)
  g_beta_etaj=g_beta_etaj_list$g
  
  
  #Score Calculation (Beta, etaj)
  #S_1 (tilt_ratio*\nabla_{\beta}log f(y_i;\bar{\beta},\bar{\gamma_j}))
  ll=function(g){Lik_tilt(g=g,x=x,y=y,model=model,tilt_ratio)}
  S= NULL
  S=-1*n*numDeriv::grad(ll, g_beta_etaj) #Likelihood is already divided by n (Likelihood multiplied by -1/n)
  S1=matrix(S[beta_indices],ncol=1)
  
  #S_2 (tilt ratio*\tilde{H}_{\beta\gamma}^{(1,j)}*(\tilde{H}_{\gamma\gamma}^{(1,j)})^{-1}\nabla_{\gamma}log(f(y_i;\beta,\bar{\gamma_j})))
  Ht=tilt_I(betabar=betabar,etaj=etaj,eta1=eta1,x=x,y=y,eta_indices=eta_indices,model=model) #Already multiplied by -1/n
  Ht_ee = Ht[eta_indices,eta_indices]
  Ht_be = Ht[beta_indices,eta_indices]
  S2=Ht_be%*%solve(Ht_ee)%*%matrix(S[eta_indices],ncol=1)
  #print("Ht_ee")
  #print(Ht_ee)
  #print("Ht_be")
  #print(Ht_be)
  #print("solve(Ht_ee)")
  #print(solve(Ht_ee))
  
  S=S1-S2
  S
}

#x and y within site
#' @keywords internal
check_U = function(beta, betabar, eta_mat, site, x,y,k,eta_indices,model){#Calculates the sum over the tilted score functions across all sites j=1,...,K
  #print(beta)
  #k is the local site?
  nsite=table(site)
  n=nrow(x)
  K = nrow(eta_mat)
  l_beta= length(beta)
  l_eta = ncol(eta_mat)
  l_tot= l_beta+l_eta
  temp = 0
  for(i in 1:K){
    temp = temp +tilt_es(beta=beta,betabar=betabar,etaj=eta_mat[i,],eta1=eta_mat[k,],x=x,y=y,eta_indices=eta_indices,model=model)
  }
  temp/(K*n)
}


#' @import plyr
#' @keywords internal
S= function(beta,eta_mat, x, y,eta_indices,model,site){ #Total sum over all S_j(\bar{\beta},\bar{\gamma}_j)
  K = nrow(eta_mat)
  unique_site=sort(unique(site))
  n_vec=count(site)
  n=length(site)
  temp = 0
  for(i in unique_site){
    #print(i)
    temp=temp +(regular_es(beta=beta,eta=eta_mat[i,],x=x[site==i,],y=y[site==i],eta_indices=eta_indices,model=model)*n_vec[i,2]) #Regular_es already accounts for n
  }
  temp/(n)
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
