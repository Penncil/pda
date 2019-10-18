# Rcpp::sourceCpp('engine/src/rcpp_coxph.cpp')
Rcpp::sourceCpp('../ODACO/src/rcpp_coxph.cpp')

#' @useDynLib ODACO
#' @import stats
#' @import Rcpp
#' @import survival 
#' @title ODACO: One-shot Distributed Algorithm of COx regression
#' 
#' @description  Fit a  Distributed Cox regression that only requires one-shot communication between 
#'               the local and remote sites, by optimizing Jordan's surrogate log-likelihood (Jordan2018JASA) function.
#' @usage DistCox(mydata, id.local=1, init_est = 'local',  output.ODACO1=F, strat=F, control=list(maxit=100), verbose = F)
#' @author Chongliang Luo, Rui Duan, Yong Chen
#' 
#' @param mydata Data frame of the data from all sites, with the first three columns as site id, event time, 
#'                 event indicator (TRUE=event observed, FALSE=censored) and the rest columns as covariates.
#' @param init_est Character or numeric vector, default 'local' which uses local est as initial beta_bar.
#' @param output.ODACO1 Logical - output ODACO1, which the surrogate L is 1st order approximation? Default FALSE.
#' @param strat   Logical, TRUE=assume heterogeneous baseline hazard in each site, FALSE=the same. Default FALSE. 
#' @param control control options passed to \code{optim()}, e.g. list(maxit=100).
#' @param verbose Logical - show running details?
#'
#' @return List with component:
#' @return \item{beta_bar}{init coef est}
#' @return \item{beta_tilde}{DistCox coef est by optimizing the surrogate L by 2nd order approximation}
#' @return \item{var_tilde}{var of beta_tilde calculated from hessian}
#' @return \item{grad_LN_L1_beta_bar}{gradient from remote machines}
#' @return \item{hess_LN_L1_beta_bar}{hessian  from remote machines}
#' @return \item{beta_tilde1}{DistCox coef est by optimizing the surrogate L by 1st order approximation}
#' @return \item{var_tilde1}{var of beta_tilde1 calculated from hessian}
##' @return \item{sol}{an \code{optim()} output, containing beta_tilde}
##' @return  \item{sol1}{an \code{optim()} output, containing beta_tilde1}
#' @details The surrogate log-likelihood is the logL of local machine plus 
#'          the 1st and 2nd order approximation of the logL of remote machines
#'          Ltilde(beta) = L1(beta) + sum(beta * grad(LN(beta_bar) - L1(beta_bar)))+ t(beta) %*% hessian(LN(beta_bar) - L1(beta_bar)) %*% beta / 2.
#'          It's optimized by \code{optim()}. The log-L, gradients and hessian are written in Rcpp.
#'          
#'          
#' @references Jordan, Michael I., Jason D. Lee, and Yun Yang. "Communication-efficient distributed statistical inference." JASA (2018).
#' @examples
#' library(survival)
#' data(lung, package = 'survival')
#' table(lung$inst)
#' # set inst#1 as the loal site
#' # choose age and sex as covariates
#' lung1 <- lung[, c('inst', 'time', 'status', 'age', 'sex')]
#' lung1$status <- lung1$status==2
#' 
#' # proposed DistCox, use likelihood from local and gradient and hessian from remote
#' fit.ODAC <- DistCox(mydata = lung1,
#'                         id.local = 1,
#'                         init_est = 'local',
#'                         strat = F,
#'                         output.ODACO1 = T)
#' # use pooled data, survival::coxph()  
#' fit.pool <- coxph(Surv(time, status) ~ age+sex, data=lung1)
#' # compare: DistCox (beta_tilde and beta_tilde1) obtains coef est better than using local only, closer to use pooled data
#' rbind(fit.ODAC$beta_bar,
#'       fit.ODAC$beta_tilde,
#'       fit.ODAC$beta_tilde1, 
#'       fit.pool$coef) 
#' # compare s.d. estimates
#' rbind(sqrt(fit.ODAC$var_tilde),
#'       sqrt(fit.ODAC$var_tilde1),
#'       summary(fit.pool)$coef[,3] )
#'       
#' ######## assume heterogeneous baseline hazard  ##########
#' ## set 3 fake sites       
#'  lung2=lung[,1:5]
#' lung2$inst[lung$inst>=2] = 2
#' lung2$inst[lung$inst>=12] = 3
#' lung2$inst[is.na(lung2$inst)] <- 3
#' table(lung2$inst)
#' lung2$status <- lung2$status == 2
#' fit.ODACH.3 = DistCox(mydata=lung2, id.local = 3, init_est = 'local', output.ODACO1 = T, verbose = T, strat = T)
#' fit.ODACH.2 = DistCox(mydata=lung2, id.local = 2, init_est = 'local', output.ODACO1 = T, verbose = T, strat = T)
#' fit.ODACH.1 = DistCox(mydata=lung2, id.local = 1, init_est = 'local', output.ODACO1 = T, verbose = T, strat = T)
#' fit.strCox <- coxph(Surv(time, status) ~ age+sex +strata(inst), data=lung2 )
#' rbind(fit.ODACH.1$beta_bar, 
#'       fit.ODACH.2$beta_bar, 
#'       fit.ODACH.3$beta_bar, 
#'       fit.ODACH.1$beta_tilde, 
#'       fit.ODACH.2$beta_tilde,
#'       fit.ODACH.3$beta_tilde
#'       fit.strCox$coef,
#'       lung_cox_pkg$coef)
#'
#' # par(mfrow=c(1,3))
#' # for(ii in 1:3){
#' #   fit.cox <- coxph(Surv(time, status) ~ age+sex, data=lung2, subset = inst==ii) 
#' #   plot(survfit(fit.cox), fun='cumhaz', main=paste0("cumulative hazard of inst ", ii), xlim=c(0,1100), ylim=c(0,3.5))
#' # }
#' 
#' @export
DistCox <- function(# local_data,  # 
                    # all_data,     
                    mydata,  # col1=id.site, col2=t_surv, col3=ind_event, col4+=X
                    id.local = 1,
                    init_est = 'local',         
                    # hasTies = NULL,            
                    verbose = F,
                    output.ODACO1 = FALSE,
                    strat = F, # T=easy ODAC, allow heterogeneous baseline haz in sites
                    control = list(maxit=100)){
  id.site <- mydata[,1]
  mydata_local <- mydata[id.site==id.local,]
  make.list <- function(mydata) list(t_surv=mydata[,2], ind_event=mydata[,3], X=as.matrix(mydata[,-c(1:3)]))
  
  local_data = all_data = list()
  local_data <- make.list(mydata_local)
  # local_data$X <- as.matrix(local_data$X)
  # all_data$X <- as.matrix(all_data$X)
  all_data <- make.list(mydata)
  
  px <- ncol(all_data$X)
  # if(is.null(hasTies)) 
  hasTies <- anyDuplicated(all_data$t_surv) > 0

  # if init_est not provided, use Cox reg of the local machine (#1) as init est beta_bar
  if(is.numeric(init_est)){
    beta_bar <- init_est
  } else {
    if(init_est!='local') warning("wrong init_est:  use local")
    sol_l <- my_coxph(local_data)
    beta_bar <- sol_l$par 
  }
    
  # logL and its gradient are written in rcpp
  if(hasTies){
    cox_fun_1_logL <- function(beta) rcpp_coxph_logL_efron(beta, time = local_data$t_surv, event = local_data$ind_event, z = local_data$X) / length(local_data$t_surv)
    cox_fun_1_logL_gradient <- function(beta) rcpp_coxph_logL_gradient_efron(beta, time = local_data$t_surv, event = local_data$ind_event, z = local_data$X) / length(local_data$t_surv)
    cox_fun_1_logL_hess <- function(beta) rcpp_coxph_logL_hessian(beta, time = local_data$t_surv, event = local_data$ind_event, z = local_data$X) / length(local_data$t_surv)
    
    if(strat==F){  # ODAC, assuming common baseline haz
      cox_fun_N_logL_gradient <- function(beta) rcpp_coxph_logL_gradient_efron(beta, time = all_data$t_surv, event = all_data$ind_event, z = all_data$X) / length(all_data$t_surv)
      grad_LN_L1_beta_bar <- cox_fun_N_logL_gradient(beta_bar) - cox_fun_1_logL_gradient(beta_bar)
      # also compute 2nd order in surrogate L
      cox_fun_N_logL_hess <- function(beta) rcpp_coxph_logL_hessian(beta, time = all_data$t_surv, event = all_data$ind_event, z = all_data$X) / length(all_data$t_surv)
      hess_LN_L1_beta_bar <- matrix(cox_fun_N_logL_hess(beta_bar) - cox_fun_1_logL_hess(beta_bar), px, px)
    } else {       # easy ODAC, each site use their own data to calculate grad, hess!
      grad_LN_L1_beta_bar <- 0
      hess_LN_L1_beta_bar <- 0
      for(ik in setdiff(id.site, id.local)){
        di <- make.list(mydata[id.site==ik,])
        grad_LN_L1_beta_bar <- grad_LN_L1_beta_bar + 
          rcpp_coxph_logL_gradient_efron(beta_bar, time = di$t_surv, event = di$ind_event, z = di$X) / length(di$t_surv)
        hess_LN_L1_beta_bar <- hess_LN_L1_beta_bar +
          rcpp_coxph_logL_hessian(beta_bar, time = di$t_surv, event = di$ind_event, z = di$X) / length(di$t_surv)
      }
      hess_LN_L1_beta_bar <- matrix(hess_LN_L1_beta_bar, px, px)
    }
  } else {
    cox_fun_1_logL <- function(beta) rcpp_coxph_logL(beta, time = local_data$t_surv, event = local_data$ind_event, z = local_data$X) / length(local_data$t_surv)
    cox_fun_1_logL_gradient <- function(beta) rcpp_coxph_logL_gradient(beta, time = local_data$t_surv, event = local_data$ind_event, z = local_data$X) / length(local_data$t_surv)
    cox_fun_1_logL_hess <- function(beta) rcpp_coxph_logL_hessian(beta, time = local_data$t_surv, event = local_data$ind_event, z = local_data$X) / length(local_data$t_surv)
    
    if(strat==F){
      cox_fun_N_logL_gradient <- function(beta) rcpp_coxph_logL_gradient(beta, time = all_data$t_surv, event = all_data$ind_event, z = all_data$X) / length(all_data$t_surv)
      grad_LN_L1_beta_bar <- cox_fun_N_logL_gradient(beta_bar) - cox_fun_1_logL_gradient(beta_bar)
      # also compute 2nd order in surrogate L
      cox_fun_N_logL_hess <- function(beta) rcpp_coxph_logL_hessian(beta, time = all_data$t_surv, event = all_data$ind_event, z = all_data$X) / length(all_data$t_surv)
      hess_LN_L1_beta_bar <- matrix(cox_fun_N_logL_hess(beta_bar) - cox_fun_1_logL_hess(beta_bar), px, px)
    } else {
      grad_LN_L1_beta_bar <- 0
      hess_LN_L1_beta_bar <- 0
      for(ik in setdiff(id.site, id.local)){
        di <- make.list(mydata[id.site==ik,])
        grad_LN_L1_beta_bar <- grad_LN_L1_beta_bar + 
          rcpp_coxph_logL_gradient(beta_bar, time = di$t_surv, event = di$ind_event, z = di$X) / length(di$t_surv)
        hess_LN_L1_beta_bar <- hess_LN_L1_beta_bar +
          rcpp_coxph_logL_hessian(beta_bar, time = di$t_surv, event = di$ind_event, z = di$X) / length(di$t_surv)
      }
      hess_LN_L1_beta_bar <- matrix(hess_LN_L1_beta_bar, px, px)
    }
  }
  

  
  # surrogate log-L and its gradient
  Ltilde1 <- function(beta) cox_fun_1_logL(beta) + sum(beta * grad_LN_L1_beta_bar)  
  Ltilde <- function(beta) cox_fun_1_logL(beta) + sum(beta * grad_LN_L1_beta_bar) + 1/2 * t(beta-beta_bar) %*% hess_LN_L1_beta_bar %*% (beta-beta_bar)
  
  Ltilde_gradient1 <- function(beta) cox_fun_1_logL_gradient(beta) + grad_LN_L1_beta_bar
  Ltilde_gradient <- function(beta) cox_fun_1_logL_gradient(beta) + grad_LN_L1_beta_bar + hess_LN_L1_beta_bar %*% (beta - beta_bar)
  
  if(output.ODACO1){
    sol1 <- optim(par = beta_bar, 
                  fn = Ltilde1,
                  gr = Ltilde_gradient1,
                  hessian = T, control = control)
  } else {
    sol1 <- list()
  }
  
  sol <- optim(par = beta_bar, 
               fn = Ltilde,
               gr = Ltilde_gradient,
               hessian = T, control = control)
  
  if(verbose) cat('conv: ODACO=', sol$convergence, 'ODACO1=', sol1$convergence, '\n')
  beta_tilde <- sol$par
  beta_tilde1 <- sol1$par
  
  # calculate var from inv hessian
  # if(output.ODACO1==TRUE){
    # if(min(svd(sol$hessian)$d)<1e-7)
      # var_tilde <- rep(NA, px)
    # else
      var_tilde <- tryCatch(diag(solve(sol$hessian))/ length(all_data$t_surv), error=function(e) NA)
  # } else{
    # var_tilde <- NA
  # }
  
  
  # if(min(svd(sol2$hessian)$d)<1e-7)   #  Error in svd(sol2$hessian) : infinite or missing values in 'x' 
    # var_tilde2 <- rep(NA, px)
  # else
    var_tilde1 <- tryCatch(diag(solve(sol1$hessian))/ length(all_data$t_surv), error=function(e) NA)
                           

  return(list(beta_bar = beta_bar,
              beta_tilde = beta_tilde,
              var_tilde  = var_tilde,
              sol = sol,
              beta_tilde1 = beta_tilde1,
              var_tilde1 = var_tilde1,
              sol1 = sol1,
              grad_LN_L1_beta_bar = grad_LN_L1_beta_bar,
              hess_LN_L1_beta_bar = hess_LN_L1_beta_bar ))
}




 
#' @title Combine and reformate data from multiple sites
#' 
#' @description  combine data (in data frame) from multiple machines together to list of time-to-event, event indicator and covariates.
#' @usage combine_data(local_data, remote_data, col_time, col_event, col_X)
#' @author Chongliang Luo, Rui Duan, Yong Chen
#' 
#' @param local_data Data frame containing time-to-event, event indicator and covariates in the local machine.
#' @param remote_data  List of data frames, each df containing time-to-event, event indicator and covariates in a remote machine.
#' @param col_time Integer, which col of the df is the time-to-event. time-to-event need to be positive, can be tied and right-censored (e.g. event=1)
#' @param col_event Integer, which col of the df is the event indicator. event=1 if event happens, event=0 if censored
#' @param col_X Vector of integers, which cols are the covariates. If not specified (i.e. col_X=0), use all col except col_time and col_event
#' 
#' @return List of local_data and all_data
#' @examples 
#' # See DistCox()
#' @export
combine_data <- function(local_data = data.frame(), 
                         remote_data = list(data.frame(), data.frame()),
                         col_time,
                         col_event,
                         col_X = 0){
  if(col_X==0) col_X <- setdiff(1:ncol(local_data), c(col_time, col_event))
  
  local_data <- list(t_surv=local_data[,col_time],   # time-to-event
                     ind_event=local_data[,col_event], # event indicator
                     X=as.matrix(local_data[,col_X]) )        # covariates
  all_data <- local_data
  for(i in 1:length(remote_data)){
    all_data$t_surv <- c(all_data$t_surv, remote_data[[i]][,col_time])
    all_data$ind_event <- c(all_data$ind_event, remote_data[[i]][,col_event])
    all_data$X <- rbind(all_data$X, as.matrix(remote_data[[i]][,col_X]))
  }
  
  return(list(local_data=local_data,
              all_data=all_data))
}
 


#' @title Fit Cox PH model 
#' 
#' @description  fit Cox reg to time-to-event data. Expected to be the same as coxph() in package 'survival'. 
#'               Optimized by optim(), with log-L and gradient written in Rcpp.
#' @usage my_coxph(data, beta_ini, hasTies = NULL)
#' @author Chongliang Luo, Rui Duan, Yong Chen
#' 
#' @param data List with component time (time-to-event), event (event indicator) and X (covariates).
#' @param beta_ini  Vector of coef for initializing \code{optim()}, if not specified, use 0.
#' @param hasTies Logical - is the time to event tied? If yes, use Efron's approximation for tie correction in the partial likelihood.
#' 
#' @return An \code{optim()} output.
#' @examples 
#' # see the example in DistCox()
#' @export
my_coxph <- function(data = list(),
                     beta_ini = c(), 
                     hasTies = NULL ){
  if(is.null(beta_ini)) beta_ini <- rep(0, ncol(data$X))
  data$X <- as.matrix(data$X)
  if(is.null(hasTies)) hasTies <- anyDuplicated(data$t_surv) > 0
  
  if(hasTies){
    cox_fun_logL <- function(beta) rcpp_coxph_logL_efron(beta, time = data$t_surv, event = data$ind_event, z = data$X) / length(data$t_surv)
    cox_fun_logL_gradient <- function(beta) rcpp_coxph_logL_gradient_efron(beta, time = data$t_surv, event = data$ind_event, z = data$X) / length(data$t_surv)
  } else {
    cox_fun_logL <- function(beta) rcpp_coxph_logL(beta, time = data$t_surv, event = data$ind_event, z = data$X) / length(data$t_surv)
    cox_fun_logL_gradient <- function(beta) rcpp_coxph_logL_gradient(beta, time = data$t_surv, event = data$ind_event, z = data$X) / length(data$t_surv)
  }
  
  sol <- optim(par = beta_ini,
                fn = cox_fun_logL,
                gr = cox_fun_logL_gradient,  
                hessian = T)

  return(sol)
}



 
