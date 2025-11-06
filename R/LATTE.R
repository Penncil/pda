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
library(glmnet)
# Set in pda()
LATTE.steps <- c('initialize', 'estimate')
LATTE.family <- 'binomial'

#' @useDynLib pda
#' @title LATTE initialize
#' @description Initialize step for LATTE (Lossless Aggregation for Treatment effect estimation)
#' 
#' @usage LATTE.initialize(ipdata, control, config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#' 
#' @return init object containing prepared data and PS model
#' @keywords internal
LATTE.initialize <- function(ipdata, control, config) {
  # Extract variables for PS model
  xvars <- colnames(ipdata)[!grepl("^outcome_", colnames(ipdata)) & 
                           !colnames(ipdata) %in% c("ID", "treatment", "index_date","site", "group")]
  xvars <- xvars[colSums(ipdata[xvars]) > 30]
  
  # Create PS model formula
  ps_formula <- as.formula(paste("treatment ~", paste(xvars, collapse = "+")))
  
  # Prepare data for PS model
  mydata_ps <- ipdata[, colnames(ipdata) %in% c(xvars, "treatment", control$outcome)]
  # Xmat <- model.matrix(ps_formula, data = mydata_ps)[,-1] # Remove intercept
  Xmat <- grab_design_matrix(data = mydata_ps, rhs_formula = ps_formula)
  Y <- mydata_ps$treatment
  
  # Fit PS model using cv.glmnet
  nfolds <- 10
  foldid <- sample(rep(seq(nfolds), length.out = length(Y)))
  set.seed(42)
  ps_fit_cv <- cv.glmnet(Xmat, Y, alpha = 1, family = "binomial", 
                        nfolds = nfolds, foldid = foldid)
  ps_fit <- glmnet(Xmat, Y, alpha = 1, family = "binomial", 
                   lambda = ps_fit_cv$lambda.min)
  
  # Calculate propensity scores
  propensityScore <- predict(ps_fit, Xmat, type = "response")
  mydata_ps$propensityScore <- propensityScore
  
  # Find optimal number of strata
  best_strata <- optimize_strata(mydata_ps, xvars)
  stratifiedPop <- get_stratified_pop(mydata_ps, nstrata = best_strata$n_strata)
  ADdata <- create_2x2_tables(
    stratifiedPop$treatment, 
    stratifiedPop[[control$outcome]], 
    stratifiedPop$stratumId
  )
  # Return initialization results
  return(list(
    prepared_data = ADdata
    # ps_model = ps_fit,
    # xvars = xvars
  ))
}

#' @useDynLib pda
#' @title LATTE LATTE.estimate
#' @description Analysis step for LATTE
#' 
#' @usage LATTE.estimate(init_data, control, config)
#' @param init_data initialization data from LATTE.initialize
#' @param control pda control data
#' @param config local site configuration
#' 
#' @return analysis results
#' @keywords internal
LATTE.estimate<- function(init_data, control, config) { 
  K <- length(control$sites)
  ADdata <- list()
  for (site_i in control$sites) {
    init_i <- pdaGet(paste0(site_i, "_initialize"), config) 
    for(j in 1:length(init_i$prepared_data)){
          ADdata[[length(ADdata) + 1]] <- init_i$prepared_data[[j]]
    }
  }

  # Fit conditional logistic regression
  results <- optimize_conditional_logistic_2x2(ADdata)
   
  # Return results
  return(list(
    coefficients = results$coefficients,
    se = results$se,
    odds_ratio = results$odds_ratios,
    ci_lower = results$ci_lower,
    ci_upper = results$ci_upper,
    convergence = results$convergence
  ))
}

# Helper Functions

#' @keywords internal
optimize_strata <- function(data, xvars, min_strata = 2, max_strata = 6) {
  best_n_strata <- 0
  best_after <- 1000
  best_smd_res <- NULL
  
  for (nstrata in min_strata:max_strata) {
    stratifiedPop <- get_stratified_pop(data, nstrata = nstrata)
    smd_res <- get_SMD(stratifiedPop = stratifiedPop, xvars = xvars)
    
    before <- sum(abs(smd_res$smd_before$SMD) > 0.2)
    after <- sum(abs(smd_res$smd_after$SMD) > 0.2)
    
    if (after < best_after) {
      best_after <- after
      best_n_strata <- nstrata
      best_smd_res <- smd_res
    }
  }
  
  return(list(
    n_strata = best_n_strata,
    smd_res = best_smd_res
  ))
}

# Include other helper functions from latte_codes.R:
# - logLik_conditional()
# - logLik_conditional_2x2()
# - optimize_conditional_logistic()
# - optimize_conditional_logistic_2x2()
# - create_2x2_tables()

get_stratified_pop = function(mydata_test, nstrata){
  rowId = c(1:length(mydata_test$treatment))
  mydata_test <- cbind(rowId, mydata_test)
  stratifiedPop <- stratifyByPs(mydata_test, numberOfStrata = nstrata)
  return (stratifiedPop)
}

stratifyByPs <- function(population, numberOfStrata = 5, stratificationColumns = c(), baseSelection = "all") {
  if (!("rowId" %in% colnames(population)))
    stop("Missing column rowId in population")
  if (!("treatment" %in% colnames(population)))
    stop("Missing column treatment in population")
  if (!("propensityScore" %in% colnames(population)))
    stop("Missing column propensityScore in population")
  ParallelLogger::logTrace("Stratifying by propensity score")
  if (nrow(population) == 0) {
    return(population)
  }
  baseSelection <- tolower(baseSelection)
  if (baseSelection == "all") {
    basePop <- population$propensityScore
  } else if (baseSelection == "target") {
    basePop <- population$propensityScore[population$treatment == 1]
  } else if (baseSelection == "comparator") {
    basePop <- population$propensityScore[population$treatment == 0]
  } else {
    stop(paste0("Unknown base selection: '", baseSelection, "'. Please choose 'all', 'target', or 'comparator'"))
  }
  if (length(basePop) == 0) {
    psStrata <- c()
  } else {
    psStrata <- unique(quantile(basePop, (1:(numberOfStrata - 1))/numberOfStrata))
  }
  attr(population, "strata") <- psStrata
  breaks <- unique(c(0, psStrata, 1))
  breaks[1] <- -1 # So 0 is included in the left-most stratum
  if (length(breaks) - 1 < numberOfStrata) {
    warning("Specified ", numberOfStrata, " strata, but only ", length(breaks) - 1, " could be created")
  }
  if (length(stratificationColumns) == 0) {
    if (length(breaks) - 1 == 1) {
      population$stratumId <- rep(1, nrow(population))
    } else {
      population$stratumId <- as.integer(as.character(cut(population$propensityScore,
                                                          breaks = breaks,
                                                          labels = 1:(length(breaks) - 1))))
    }
    return(population)
  } else {
    f <- function(subset, psStrata, numberOfStrata) {
      if (length(breaks) - 1 == 1) {
        subset$stratumId <- rep(1, nrow(subset))
      } else {
        subset$stratumId <- as.integer(as.character(cut(subset$propensityScore,
                                                        breaks = breaks,
                                                        labels = 1:(length(breaks) - 1))))
      }
      return(subset)
    }
    
    results <- plyr::dlply(.data = population,
                           .variables = stratificationColumns,
                           .drop = TRUE,
                           .fun = f,
                           psStrata = psStrata,
                           numberOfStrata = numberOfStrata)
    maxStratumId <- 0
    for (i in 1:length(results)) {
      if (nrow(results[[i]]) > 0) {
        if (maxStratumId != 0)
          results[[i]]$stratumId <- results[[i]]$stratumId + maxStratumId + 1
        maxStratumId <- max(results[[i]]$stratumId)
      }
    }
    result <- do.call(rbind, results)
    return(result)
  }
}
get_SMD <-function(stratifiedPop, xvars){
  smd_before <- GetSMD(stratifiedPop[, xvars], stratifiedPop$treatment)
  weight_mat=Compute_weight(stratifiedPop)
  stratifiedPop=stratifiedPop%>%left_join(weight_mat,by="rowId")
  smd_after <- GetSMD(stratifiedPop[, xvars], stratifiedPop$treatment,stratifiedPop$weight)
  
  return(list(smd_before=smd_before,smd_after=smd_after))
}
GetSMD <- function(data, treat, weights = NULL, std = TRUE){
  table <- col_w_smd(data, treat, weights, std)
  smd = t(data.frame(as.list(table)))
  smd = data.frame(row.names(smd), smd)
  rownames(smd) <- NULL
  colnames(smd) = c("Xvars", "SMD")
  
  return(smd)
}


Compute_weight <- function(data) {
  stratumSize <- data %>%
    group_by(stratumId, treatment) %>%
    count() %>%
    ungroup()

  w <- stratumSize %>%
    mutate(weight = 1 / n) %>%
    inner_join(data, by = c("stratumId", "treatment"), multiple = "all") %>%
    select(rowId, treatment, weight)

  wSum <- w %>%
    group_by(treatment) %>%
    summarize(wSum = sum(weight, na.rm = TRUE)) %>%
    ungroup()

  w_final <- w %>%
    inner_join(wSum, by = "treatment") %>%
    mutate(weight = weight / wSum) %>%
    select(rowId, treatment, weight)

  return(w_final[, c(1, 3)])
}


optimize_conditional_logistic_2x2 <- function(tables, start_beta = 0.1) {  # Changed default start_beta
  result <- optim(par = start_beta, 
                 fn = logLik_conditional_2x2, 
                 tables = tables, 
                 method = "BFGS",
                 control = list(
                   maxit = 1000,
                   reltol = 1e-8,
                   ndeps = 1e-6
                 ),
                 hessian = TRUE)
  
  if (result$convergence != 0) {
    warning("Optimization did not converge. Results may be unreliable.")
  }
  
  # Check if hessian is invertible
  hessian_inv <- try(solve(result$hessian), silent = TRUE)
  if (inherits(hessian_inv, "try-error")) {
    warning("Hessian matrix is not invertible. Standard errors may be unreliable.")
    se <- NA
  } else {
    se <- sqrt(diag(hessian_inv))
  }
  
  odds_ratios <- exp(result$par)
  ci_lower <- exp(result$par - 1.96 * se)
  ci_upper <- exp(result$par + 1.96 * se)
  
  list(
    coefficients = result$par,
    se = se,
    odds_ratios = odds_ratios,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    log_likelihood = -result$value,
    convergence = result$convergence,
    message = result$message
  )
}

# Optimize conditional logistic regression
logLik_conditional_2x2 <- function(beta, tables) {
  loglik <- 0

  for (table in tables) {
    # Using standard notation for 2x2 tables
    a <- table[1,1]  # exposed cases
    b <- table[1,2]  # unexposed cases
    c <- table[2,1]  # exposed controls
    d <- table[2,2]  # unexposed controls

    n1 <- a + b      # total cases
    n2 <- c + d      # total controls
    m1 <- a + c      # total exposed
    
    # Skip uninformative strata
    # if (n1 == 0 || n2 == 0 || m1 == 0 || m1 == (n1 + n2)) {
    #   next
    # }
    if (n1 == 0 || n2 == 0) {  
    # Instead of skipping, assign a small likelihood contribution  
        loglik <- loglik + log(1e-10)  # Avoid -Inf
        next
    }


    # Contribution to the log-likelihood
    numerator <- beta * a
    
    # Calculate denominator using log-sum-exp for numerical stability
    lower_k <- max(0, n1 - (n1 + n2 - m1))
    upper_k <- min(n1, m1)
    
    log_terms <- sapply(lower_k:upper_k, function(k) {
      lgamma(m1 + 1) - lgamma(k + 1) - lgamma(m1 - k + 1) + 
      lgamma(n1 + n2 - m1 + 1) - lgamma(n1 - k + 1) - lgamma(n2 - (m1 - k) + 1) + 
      beta * k
    })
    
    max_log_term <- max(log_terms)
    denominator <- max_log_term + log(sum(exp(log_terms - max_log_term)))
    
    loglik <- loglik + (numerator - denominator)
  }
  print(loglik)
  return(-loglik)  # Return negative for minimization
}
optimize_conditional_logistic_2x2 <- function(tables, start_beta = 0.1) {  # Changed default start_beta
  result <- optim(par = start_beta, 
                 fn = logLik_conditional_2x2, 
                 tables = tables, 
                 method = "BFGS",
                 control = list(
                   maxit = 1000,
                   reltol = 1e-8,
                   ndeps = 1e-6
                 ),
                 hessian = TRUE)
  
  if (result$convergence != 0) {
    warning("Optimization did not converge. Results may be unreliable.")
  }
  
  # Check if hessian is invertible
  hessian_inv <- try(solve(result$hessian), silent = TRUE)
  if (inherits(hessian_inv, "try-error")) {
    warning("Hessian matrix is not invertible. Standard errors may be unreliable.")
    se <- NA
  } else {
    se <- sqrt(diag(hessian_inv))
  }
  
  odds_ratios <- exp(result$par)
  ci_lower <- exp(result$par - 1.96 * se)
  ci_upper <- exp(result$par + 1.96 * se)
  
  list(
    coefficients = result$par,
    se = se,
    odds_ratios = odds_ratios,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    log_likelihood = -result$value,
    convergence = result$convergence,
    message = result$message
  )
}
