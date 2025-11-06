
logLik_conditional <- function(beta, X, Y, strata) {
  n_strata <- length(unique(strata))
  loglik <- 0
  
  # Add bounds checking for beta
  if (any(abs(beta) > 50)) {
    return(.Machine$double.xmax)
  }
  
  for (s in 1:n_strata) {
    stratum_indices <- which(strata == s)
    X_s <- X[stratum_indices, drop = FALSE]
    Y_s <- Y[stratum_indices]
    
    n_cases <- sum(Y_s)
    if (n_cases == 0 || n_cases == length(Y_s)) {
      next
    }
    
    X_sum <- sum(X_s[Y_s == 1])
    numerator <- beta * X_sum
    
    n_total <- length(Y_s)
    X_total <- sum(X_s)
    
    # Early return for numerical instability
    if (is.nan(X_sum) || is.infinite(X_sum)) {
      return(.Machine$double.xmax)
    }
    
    # Modified log-sum-exp with additional checks
    k_values <- 0:min(n_cases, X_total)
    log_terms <- numeric(length(k_values))
    
    for (i in seq_along(k_values)) {
      k <- k_values[i]
      # Use log1p for better numerical stability when needed
      log_choose1 <- lchoose(X_total, k)
      log_choose2 <- lchoose(n_total - X_total, n_cases - k)
      
      if (is.infinite(log_choose1) || is.infinite(log_choose2)) {
        next
      }
      
      log_terms[i] <- log_choose1 + log_choose2 + beta * k
    }
    
    # Remove any infinite values
    log_terms <- log_terms[is.finite(log_terms)]
    
    if (length(log_terms) == 0) {
      return(.Machine$double.xmax)
    }
    
    max_log_term <- max(log_terms)
    denominator <- max_log_term + log(sum(exp(log_terms - max_log_term)))
    
    contribution <- numerator - denominator
    if (is.finite(contribution)) {
      loglik <- loglik + contribution
    }
  }
  
  if (is.nan(loglik) || is.infinite(loglik)) {
    return(.Machine$double.xmax)
  }
  
  return(-loglik)
}

optimize_conditional_logistic <- function(X, Y, strata, start_beta = NULL) {
  if (is.null(start_beta)) {
    start_beta <- rep(0.1, ncol(X))  # Changed from -1 to 0.1 for better initial values
  }
  
  # Add bounds to optimization
  result <- optim(par = start_beta, 
                 fn = logLik_conditional, 
                 X = X, Y = Y, strata = strata,
                 method = "BFGS",
                 control = list(
                   maxit = 1000,
                   reltol = 1e-8,
                   ndeps = rep(1e-6, length(start_beta))
                 ),
                 hessian = TRUE)
  
  if (result$convergence != 0) {
    warning("Optimization did not converge. Results may be unreliable.")
  }
  
  # Check if hessian is invertible
  hessian_inv <- try(solve(result$hessian), silent = TRUE)
  if (inherits(hessian_inv, "try-error")) {
    warning("Hessian matrix is not invertible. Standard errors may be unreliable.")
    se <- rep(NA, length(result$par))
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



create_2x2_tables <- function(X, Y, strata) {
  unique_strata <- unique(strata)
  tables <- list()

  for (s in unique_strata) {
    # Filter data for current stratum
    X_s <- X[strata == s]
    Y_s <- Y[strata == s]

    # Create 2x2 table
    print(X)
    print(Y)
    print(strata)
    table <- matrix(0, nrow = 2, ncol = 2)
    table[1, 1] <- sum(X_s == 1 & Y_s == 1) # a: Exposed cases
    table[1, 2] <- sum(X_s == 0 & Y_s == 1) # b: Exposed controls
    table[2, 1] <- sum(X_s == 1 & Y_s == 0) # c: Unexposed cases
    table[2, 2] <- sum(X_s == 0 & Y_s == 0) # d: Unexposed controls

    # Add table to list
    tables[[as.character(s)]] <- table
  }

  return(tables)
}

# create_2x2_tables <- function(treatment, outcome, stratum) {
#   # Create a data frame with the necessary columns
#   data <- data.frame(treatment = treatment, outcome = outcome, stratum = stratum)
  
#   # Get unique strata
#   unique_strata <- unique(stratum)
  
#   # Initialize list to store tables
#   tables <- list()
  
#   # Create a 2x2 table for each stratum
#   for (s in unique_strata) {
#     # Subset data for this stratum
#     stratum_data <- data[data$stratum == s, ]
    
#     # Skip if there's no variation in treatment or outcome
#     if (length(unique(stratum_data$treatment)) < 2 || length(unique(stratum_data$outcome)) < 2) {
#       next
#     }
    
#     # Create the 2x2 table
#     table_2x2 <- matrix(0, nrow = 2, ncol = 2)
    
#     # Fill in the cells
#     # [exposed cases, unexposed cases]
#     # [exposed controls, unexposed controls]
#     table_2x2[1,1] <- sum(stratum_data$treatment == 1 & stratum_data$outcome == 1)
#     table_2x2[1,2] <- sum(stratum_data$treatment == 0 & stratum_data$outcome == 1)
#     table_2x2[2,1] <- sum(stratum_data$treatment == 1 & stratum_data$outcome == 0)
#     table_2x2[2,2] <- sum(stratum_data$treatment == 0 & stratum_data$outcome == 0)
    
#     # Add to list if the table is informative
#     tables[[as.character(s)]] <- table
#   }
  
#   return(tables)
# }

