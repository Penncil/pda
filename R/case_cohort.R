
## prepare calculation for case-cohort  design
prepare_case_cohort <- function(ipdata, method, full_cohort_size){
  # for each site, pre-calculate the failure time points, the risk sets, and the respective weights

  covariate <- as.matrix(ipdata[, -c(1:3)])
  # find over which position lies the failure times
  failure_position <- which(ipdata$status == 1)
  # find failure times
  failure_times <- ipdata$time[which(ipdata$status == 1)]
  # the number of failures
  failure_num <- length(failure_times)
  # full_cohort_size = nrow(ipdata) # ?
  Barlow_IPW <- full_cohort_size / sum(ipdata$subcohort)
 
  risk_size <- 0
  risk_sets <- as.list(rep(NA, failure_num ))
  risk_set_weights <- as.list(rep(NA, failure_num ))
  for(j in 1:failure_num){
    my_risk_set1 <- which((ipdata$subcohort == 1) & (ipdata$time >= failure_times[j]))
    risk_size <- risk_size + length(my_risk_set1)
    if(method == "Prentice"){
      my_weight1 <- rep(1, length(my_risk_set1))
    }
    if(method == "Barlow"){
      my_weight1 <- rep(Barlow_IPW, length(my_risk_set1))
    }
    if(ipdata$subcohort[which(ipdata$time == failure_times[j])] == 0){
      my_risk_set2 <- which(ipdata$time == failure_times[j])
      my_weight2 <- 1
    }else{
      my_risk_set2 <- c()
      my_weight2 <- c()
      if(method == "Barlow"){
        my_weight1[which(ipdata$time[my_risk_set1] == failure_times[j])] <- 1 
      }
    }
    risk_sets[[j]] <- c(my_risk_set1, my_risk_set2)
    risk_set_weights[[j]] <- c(my_weight1, my_weight2)
  } 
  
  return(list(# ipdata = ipdata,
              full_cohort_size = full_cohort_size,
              covariate = covariate,
              failure_position = failure_position,
              failure_num = failure_num,
              risk_sets = risk_sets,
              risk_set_weights = risk_set_weights ))
}

# this function calculate the log pseudo-likelihood for ONE site
# cc_prep is the output of prepare_case_cohort()
log_plk <- function(beta, cc_prep){
  X = cc_prep$covariate # ipdata[,-c(1:2)]
  failure_position = cc_prep$failure_position 
  failure_num = cc_prep$failure_num 
  risk_sets = cc_prep$risk_sets 
  risk_set_weights = cc_prep$risk_set_weights
  
  numerator_terms <- c(X[failure_position, ] %*% beta) 
  res <- sum(numerator_terms)
  for(j in 1:failure_num){
    temp_term <- sum(c(exp(X[risk_sets[[j]], ] %*% beta )) * risk_set_weights[[j]])
    if(temp_term > 0){
      res <- res - log(temp_term)
    }else{
      res <- res - numerator_terms[j]
    }
  }
  return(res)
}



# this function calculate the gradient of log pseudo-likelihood for ONE site
# cc_prep is the output of prepare_case_cohort()
grad_plk <- function(beta, cc_prep){
  X = cc_prep$covariate # ipdata[,-c(1:2)]
  failure_position = cc_prep$failure_position 
  failure_num = cc_prep$failure_num 
  risk_sets = cc_prep$risk_sets 
  risk_set_weights = cc_prep$risk_set_weights
  
  res = colSums(X[failure_position,])
  for(j in 1:failure_num){
    if(length(risk_sets[[j]]) > 1){
      temp_scalar <- c(exp(X[risk_sets[[j]], ] %*% beta )) * risk_set_weights[[j]]
      res <- res - apply(sweep(X[risk_sets[[j]], ], 1, temp_scalar, "*"), 2, sum) / sum(temp_scalar) 
    }else{
      res <- res - X[risk_sets[[j]], ] * risk_set_weights[[j]]
    }
  }
  return(res)
}


# this function calculate the Hessian of log pseudo-likelihood for ONE site
# cc_prep is the output of prepare_case_cohort
hess_plk <- function(beta, cc_prep){
  X = cc_prep$covariate # ipdata[,-c(1:2)]
  failure_position = cc_prep$failure_position 
  failure_num = cc_prep$failure_num 
  risk_sets = cc_prep$risk_sets 
  risk_set_weights = cc_prep$risk_set_weights
  
  d <- length(beta)
  if(length(risk_sets[[1]]) > 1){
    temp_scalar <- c(exp(X[risk_sets[[1]], ] %*% beta )) * risk_set_weights[[1]]
    temp_vec <- apply(sweep(X[risk_sets[[1]], ], 1, temp_scalar, "*"), 2, sum)
    temp_mat <- sweep(X[risk_sets[[1]], ], 1, sqrt(temp_scalar), "*")
    res <- temp_vec %*% t(temp_vec)  / (sum(temp_scalar))^2  - crossprod(temp_mat) / sum(temp_scalar)
  }else{
    res <- matrix(0, d, d)
  }
  if(failure_num > 1){
    for(j in 2:failure_num){
      if(length(risk_sets[[j]]) > 1){
        temp_scalar <- c(exp(X[risk_sets[[j]], ] %*% beta )) * risk_set_weights[[j]]
        temp_vec <- apply(sweep(X[risk_sets[[j]], ], 1, temp_scalar, "*"), 2, sum)
        temp_mat <- sweep(X[risk_sets[[j]], ], 1, sqrt(temp_scalar), "*")
        res <- res - crossprod(temp_mat) / sum(temp_scalar) + temp_vec %*% t(temp_vec)  / (sum(temp_scalar))^2
      }
    }
  }
  
  return(res)
}

