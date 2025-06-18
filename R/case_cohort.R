
## prepare calculation for case-cohort design at ONE site
# the purpose of this function is to "pre-calculate" the weight before calculating the log-likelihood
# this would accelerate the subsequent calculation of log-likelihood
# currently, we only provide Prentice weight; more options will be provided later
## this is Yudong's weight_CC() in functions_CC_1.R, can take multi-site data, or single-site as a list of length 1
## data_list contains list of ipdata, with columns: time, status, subcohort, and covariates
prepare_case_cohort <- function(data_list, method, full_cohort_size){
  # for each site, pre-calculate the failure time points, the risk sets, and the respective weights
  # also, remove those sites with zero events
  site_to_remove <- c()
  K <- length(full_cohort_size)
  failure_num <- rep(NA, K)
  failure_times <- as.list(rep(NA, K))
  risk_sets <- as.list(rep(NA, K))
  risk_set_weights <- as.list(rep(NA, K))
  covariate_list <- as.list(rep(NA, K))
  failure_position <- as.list(rep(NA, K))
  for(k in 1:K){
    # prepare a list for covariates in matrix format so as to speed up computation of log partial likelihood, gradient, and hessian
    covariate_list[[k]] <- as.matrix(data_list[[k]][, -c(1:3)]) 
    # find over which position lies the failure times
    failure_position[[k]] <- which(data_list[[k]]$status == 1)
    # find failure times
    failure_times[[k]] <- data_list[[k]]$time[which(data_list[[k]]$status == 1)]
    # the number of failures
    failure_num[k] <- length(failure_times[[k]])
    
    if(failure_num[k] == 0){
      site_to_remove <- c(site_to_remove, k)
    }else{
      risk_size <- 0
      temp_risk <- as.list(rep(NA, failure_num[k]))
      temp_weight <- as.list(rep(NA, failure_num[k]))
      for(j in 1:failure_num[k]){
        my_risk_set1 <- which((data_list[[k]]$subcohort == 1) & (data_list[[k]]$time >= failure_times[[k]][j]))
        risk_size <- risk_size + length(my_risk_set1)
        if(method == "Prentice"){
          my_weight1 <- rep(1, length(my_risk_set1))
          if(data_list[[k]]$subcohort[which(data_list[[k]]$time == failure_times[[k]][j])] == 0){
            my_risk_set2 <- which(data_list[[k]]$time == failure_times[[k]][j])
            my_weight2 <- 1
          }else{
            my_risk_set2 <- c()
            my_weight2 <- c()
          }
        } # else if (method == "Barlow")
        temp_risk[[j]] <- c(my_risk_set1, my_risk_set2)
        temp_weight[[j]] <- c(my_weight1, my_weight2)
      }
      risk_sets[[k]] <- temp_risk
      risk_set_weights[[k]] <- temp_weight
      if(risk_size == 0){
        site_to_remove <- c(site_to_remove, k)
      }
    }
  }
  
  if(length(site_to_remove) > 0){
    data_list <- data_list[-site_to_remove]
    full_cohort_size <- full_cohort_size[-site_to_remove]
    failure_num <- failure_num[-site_to_remove]
    failure_times <- failure_times[-site_to_remove]
    failure_position <- failure_position[-site_to_remove]
    risk_sets <- risk_sets[-site_to_remove]
    risk_set_weights <- risk_set_weights[-site_to_remove]
    covariate_list <- covariate_list[-site_to_remove]
    K <- K - length(site_to_remove)
  }
  
  return(list(# data_list = data_list,
              full_cohort_size = full_cohort_size,
              covariate_list = covariate_list,
              failure_position = failure_position,
              failure_num = failure_num,
              risk_sets = risk_sets,
              risk_set_weights = risk_set_weights,
              site_num=K))
}

## below only take input single-site ipdata...
# prepare_case_cohort <- function(ipdata, full_cohort_size, method){
#   # for each site, pre-calculate the failure time points, the risk sets, and the respective weights
#   # also, remove those sites with zero events
#   # for each site, pre-calculate the failure time points, the risk sets, and the respective weights
#   # ipdata columns: time, status, subcohort, and covariates
#   covariate <- as.matrix(ipdata[, -c(1:3)])
#   # find over which position lies the failure times
#   failure_position <- which(ipdata$status == 1)
#   # find failure times
#   failure_times <- ipdata$time[which(ipdata$status == 1)]
#   # the number of failures
#   failure_num <- length(failure_times) 
#   
#   risk_sets <- as.list(rep(NA, failure_num))
#   risk_set_weights <- as.list(rep(NA, failure_num))
#   
#   # if have any events
#   if(failure_num > 0) { 
#     for(j in 1:failure_num){
#       my_risk_set1 <- which((ipdata$subcohort == 1) & (ipdata$time >= failure_times[j]))
#       # risk_size <- risk_size + length(my_risk_set1)
#       if(method == "Prentice"){
#         my_weight1 <- rep(1, length(my_risk_set1))
#         if(ipdata$subcohort[which(ipdata$time == failure_times[j])] == 0){
#           my_risk_set2 <- which(ipdata$time == failure_times[j])
#           my_weight2 <- 1
#         }else{
#           my_risk_set2 <- c()
#           my_weight2 <- c()
#         }
#       } # else if(method == "Barlow")
#       risk_sets[[j]] <- c(my_risk_set1, my_risk_set2)
#       risk_set_weights[[j]] <- c(my_weight1, my_weight2)
#     }
#   }
#   
#   return(list(full_cohort_size = full_cohort_size,
#               covariate = covariate,
#               failure_position = failure_position,
#               failure_num = failure_num,
#               risk_sets = risk_sets,
#               risk_set_weights = risk_set_weights ))
# }
 

# this function calculate the log pseudo-likelihood for ONE site
# cc_prep is the output of prepare_case_cohort()
log_plk <- function(beta, cc_prep, site_num) {
  eta <- cc_prep$covariate_list[[site_num]] %*% beta
  exp_eta <- exp(eta)
  res <- sum(eta[cc_prep$failure_position[[site_num]]])
  
  for (j in 1:cc_prep$failure_num[site_num]) {
    idx <- cc_prep$risk_sets[[site_num]][[j]]
    weights <- cc_prep$risk_set_weights[[site_num]][[j]]
    res <- res - log(sum(exp_eta[idx] * weights) + 1e-12)
  }
  return(res)
}

# this function calculate the gradient of log pseudo-likelihood for ONE site
# cc_prep is the output of prepare_case_cohort()
grad_plk <- function(beta, cc_prep, site_num) {
  X <- cc_prep$covariate_list[[site_num]]
  eta <- X %*% beta
  exp_eta <- exp(eta)
  
  grad <- colSums(X[cc_prep$failure_position[[site_num]], , drop = FALSE])
  
  for (j in 1:cc_prep$failure_num[site_num]) {
    idx <- cc_prep$risk_sets[[site_num]][[j]]
    weights <- cc_prep$risk_set_weights[[site_num]][[j]]
    temp_w <- exp_eta[idx] * weights
    denom <- sum(temp_w)
    weighted_X <- sweep(X[idx, , drop = FALSE], 1, temp_w, '*')
    grad <- grad - colSums(weighted_X) / denom
  }
  return(grad)
}



# this function calculate the Hessian of log pseudo-likelihood for ONE site
# cc_prep is the output of prepare_case_cohort()
hess_plk <- function(beta, cc_prep, site_num) {
  X <- cc_prep$covariate_list[[site_num]]
  eta <- X %*% beta
  exp_eta <- exp(eta)
  d <- ncol(X)
  H <- matrix(0, d, d)
  
  for (j in 1:cc_prep$failure_num[site_num]) {
    idx <- cc_prep$risk_sets[[site_num]][[j]]
    weights <- cc_prep$risk_set_weights[[site_num]][[j]]
    temp_w <- exp_eta[idx] * weights
    denom <- sum(temp_w)
    
    X_sub <- X[idx, , drop = FALSE]
    weighted_X <- sweep(X_sub, 1, temp_w, '*')
    mean_vec <- colSums(weighted_X)
    
    sqrt_wX <- sweep(X_sub, 1, sqrt(temp_w), '*')
    
    H <- H + (tcrossprod(mean_vec) / (denom^2)) - (crossprod(sqrt_wX) / denom)
  }
  return(H)
}



# this function fits Cox PH to case-cohort (survival::cch) with the pooled multi-site data
# notice this assumes varying baseline hazard functions across sites
# cc_prep is the output of prepare_case_cohort()
#' @export
cch_pooled <- function(formula, data, subcoh='subcohort', site='site', variables_lev,
                       full_cohort_size, method = "Prentice", optim_method = "BFGS",
                       var_sandwich=T){
  n = nrow(data)
  site_uniq = unique(data[,site])
  mf <- model.frame(formula, data, xlev=variables_lev)
  
  ipdata = data.table::data.table(site=data[,site],
                                  time=as.numeric(model.response(mf))[1:n],
                                  status=as.numeric(model.response(mf))[-c(1:n)],
                                  subcohort = data[,subcoh],
                                  model.matrix(formula, mf)[,-1])
  ipdata = data.table(data.frame(ipdata))
  risk_factor = colnames(ipdata)[-c(1:4)]
  
  # notice here we allow data degeneration (e.g. missing categories in some site)
  px = ncol(ipdata)-4
  initial_beta = rep(0, px)
  names(initial_beta) = names(ipdata)[-c(1:4)]
  # pool_fun <- function(beta) sum(sapply(site_uniq, function(site_id)
  #   log_plk(beta, prepare_case_cohort(ipdata[site==site_id,-'site'], method, full_cohort_size[site_id]))))
  
  data_split <- split(ipdata, by=site, keep.by=F)
  cc_prep = prepare_case_cohort(data_split, method, full_cohort_size)
  K = cc_prep$site_num
  pool_fun <- function(beta) { 
    sum(vapply(1:K, function(i) rcpp_cc_log_plk(beta, site_num = i, 
                                                covariate_list = cc_prep$covariate_list,
                                                failure_position = cc_prep$failure_position,
                                                failure_num = cc_prep$failure_num,
                                                risk_sets = cc_prep$risk_sets,
                                                risk_set_weights = cc_prep$risk_set_weights), numeric(1)))
  }
  
  result <- optim(par = initial_beta, fn = pool_fun,
                  control = list(fnscale = -1), method = optim_method, hessian = T)
  b_pooled = result$par
  
  # calculate sandwich var estimate, degenerated data columns are given 0 coefs
  if(var_sandwich==T){
    block1 <- result$hessian
    block2 <- NULL
    data_split <- split(ipdata, ipdata$site)
    
    for(i in 1:length(site_uniq)){
      site_id <- site_uniq[i]
      ipdata_i = data_split[[i]]
      col_deg = apply(ipdata_i[,-c(1:4)],2,var)==0    # degenerated X columns...
      ipdata_i = ipdata_i[,-(which(col_deg)+4),with=F]
      # use coxph(Surv(time_in, time, status)~.) to do cch...
      precision <- min(diff(sort(ipdata_i$time))) / 2 #
      ipdata_i$time_in = 0
      ipdata_i[ipdata_i$subcohort == 0, "time_in"] <- ipdata_i[ipdata_i$subcohort == 0, "time"] - precision
      
      formula_i <- as.formula(paste("Surv(time_in, time, status) ~", paste(risk_factor[!col_deg], collapse = "+"), '+ cluster(ID)'))
      cch_i <- tryCatch(coxph(formula_i, data=cbind(ID=1:nrow(ipdata_i), ipdata_i), init=b_pooled[!col_deg], iter=0), error=function(e) NULL)
      score_resid <- resid(cch_i, type = "score")  # n x p matrix
      S_i = matrix(0, px, px)   # this is the meat in sandwich var
      S_i[!col_deg, !col_deg] <- crossprod(score_resid)
      
      block2[[i]] <- S_i
    }
    
    var <- solve(block1) %*% Reduce("+", block2) %*% solve(block1)
    result$var <- var # this is the output for variance estimates
  }
  
  return(result)
}
