require(survival)
require(data.table)
require(glmnet)
require(dplyr)
require(Matrix)
require(tibble)
library(cobalt)
library(geex)
library(numDeriv)
library(EmpiricalCalibration)

################################################################################################################################
###########################     Stratification: Conditional Logistic Regression     ############################################
################################################################################################################################

#' @keywords internal
create_2x2_tables <- function(X, Y, strata) {
  unique_strata <- unique(strata)
  tables <- list()

  for (s in unique_strata) {
    # Filter data for current stratum
    X_s <- X[strata == s]
    Y_s <- Y[strata == s]

    # Create 2x2 table
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

#' @keywords internal
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
  return(-loglik)  # Return negative for minimization
}

#' @keywords internal
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

#' @keywords internal
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
    effect_size = odds_ratios,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    log_likelihood = -result$value,
    convergence = result$convergence,
    message = result$message
  )
}


################################################################################################################################
###########################     Stratification: Conditional Poisson Regression      ############################################
################################################################################################################################ 

#' @keywords internal
create_2x2_tables_poisson <- function(stratifiedPop, outcome_id, outcome_time) {
    KSiteAD_uf = list()
    for (strat in unique(stratifiedPop$stratumId)) {
        strat_data <- stratifiedPop[stratifiedPop$stratumId == strat, ]
        
        # For Poisson, we need count data and person-time
        count_var <- outcome_id
        time_var <- outcome_time
        
        # Only include strata with variation in treatment and non-zero counts
        if (length(unique(strat_data$treatment)) > 1 && sum(strat_data[[count_var]]) > 0) {
        
        # Create the 2x2 contingency table for Poisson
        # [treated_value, treated_persontime]
        # [untreated_value, untreated_persontime]
        table_2x2 <- matrix(0, nrow = 2, ncol = 2)
        
        # Treated count
        table_2x2[1,1] <- sum(strat_data[strat_data$treatment == 1, count_var])
        
        # Treated person-time
        table_2x2[1,2] <- sum(strat_data[strat_data$treatment == 1, time_var])
        
        # Untreated count
        table_2x2[2,1] <- sum(strat_data[strat_data$treatment == 0, count_var])
        
        # Untreated person-time
        table_2x2[2,2] <- sum(strat_data[strat_data$treatment == 0, time_var])
        
        # Only add tables that are informative for Poisson regression
        # Need positive person-time in both groups and at least one event
        if (table_2x2[1,2] > 0 && table_2x2[2,2] > 0 && (table_2x2[1,1] + table_2x2[2,1]) > 0) {
            # KSiteAD_uf[[length(KSiteAD_uf) + 1]] <- table_2x2
            KSiteAD_uf[[as.character(strat)]] <- table_2x2
        }
      }
    }
    
    return(KSiteAD_uf)
}

#' @keywords internal
logLik_conditional_poisson_2x2 <- function(beta, tables) {
  loglik <- 0

  for (table in tables) {
    # For Poisson 2x2 tables:
    # [treated_value, treated_persontime]
    # [untreated_value, untreated_persontime]
    y1 <- table[1,1]  # treated count
    t1 <- table[1,2]  # treated person-time
    y0 <- table[2,1]  # untreated count
    t0 <- table[2,2]  # untreated person-time
    
    # Total events in this stratum
    total_events <- y1 + y0
    
    # Skip uninformative strata
    if (total_events == 0 || t1 == 0 || t0 == 0) {
      # Instead of skipping, assign a small likelihood contribution
      loglik <- loglik + log(1e-10)  # Avoid -Inf
      next
    }
    
    # For conditional Poisson, conditioning on total events and person-time
    # the likelihood is a function of beta, y1, and the offset ratio
    # For each stratum, we sum over all possible values of y1 from 0 to total_events
    
    # Contribution to the log-likelihood: log(P(Y1=y1 | Y1+Y0=total_events))
    # The formula uses the binomial PMF with probability adjusted by exposure and rate ratio
    
    # Numerator: probability of observed configuration
    prob <- exp(beta) * t1 / (exp(beta) * t1 + t0)
    numerator <- y1 * log(prob) + y0 * log(1 - prob)
    
    # Denominator: sum over all possible configurations
    log_terms <- sapply(0:total_events, function(k) {
      # probability for this configuration
      p <- exp(beta) * t1 / (exp(beta) * t1 + t0)
      # binomial probability
      dbinom(k, total_events, p, log = TRUE)
    })
    
    max_log_term <- max(log_terms)
    denominator <- max_log_term + log(sum(exp(log_terms - max_log_term)))
    
    # Add the contribution of this stratum
    stratum_loglik <- numerator - denominator
    loglik <- loglik + stratum_loglik
  }
  
  return(-loglik)  # Return negative for minimization
}

#' @keywords internal
# Optimization function for conditional Poisson regression
optimize_conditional_poisson_2x2 <- function(tables, start_beta = NULL) {
  if (is.null(start_beta)) {
    start_beta <- rep(0, 1)  # Starting at 0 (rate ratio = 1) is often stable
  }
  
  result <- optim(par = start_beta, 
                  fn = logLik_conditional_poisson_2x2, 
                  tables = tables, 
                  method = "BFGS",
                  hessian = TRUE)
  
  if (result$convergence != 0) {
    warning("Optimization did not converge. Results may be unreliable.")
  }
  
  se <- sqrt(diag(solve(result$hessian)))
  rate_ratios <- exp(result$par)
  ci_lower <- exp(result$par - 1.96 * se)
  ci_upper <- exp(result$par + 1.96 * se)
  
  list(
    coefficients = result$par,
    se = se,
    effect_size = rate_ratios,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    log_likelihood = -result$value,
    convergence = result$convergence
  )
}



################################################################################################################################
###################################################     IPTW      ##############################################################
################################################################################################################################ 

#' @keywords internal
trimByW <- function(propensityScore, treatment, trimFraction = 0.05) {
  # 0.05 cutoff: Stürmer  T, Rothman  KJ, Avorn  J, et al.  Treatment effects in the presence of unmeasured confounding: dealing with observations in the tails of the propensity score distribution—a simulation study. Am J Epidemiol. 2010;172(7):843–854.
  cutoffTarget <- quantile(propensityScore, trimFraction) 
  cutoffUpper <- quantile(propensityScore, 1-trimFraction)
  result <- which((propensityScore >= cutoffTarget) & ( propensityScore <= cutoffUpper) )
  return(result)
}

#' @keywords internal
getAD_IPW <- function(SiteIPD, outcome_name, formula, link = "canonical", cut_off = NULL) {
  AD = list()
  Xmat <- grab_design_matrix(data = SiteIPD, rhs_formula = formula)
  Y <- SiteIPD[[outcome_name]]
  weights <- SiteIPD$weights
  Xmat.tbl <- data.frame(Xmat)
  category_combinations <- expand.grid(lapply(Xmat.tbl, unique),
    stringsAsFactors = FALSE
  ) %>% arrange_all()
  colnames(Xmat.tbl) <- colnames(category_combinations)

  if (link == "canonical") {
    Xmat.tbl <- as_tibble(Xmat.tbl)
    cols <- colnames(Xmat.tbl)
    Xmat.tbl.weights <- cbind(Xmat.tbl, weights = weights)
    Xtable_initial <- Xmat.tbl.weights %>%
      group_by_at(.vars = cols) %>%
      summarise(n = sum(weights), .groups = "drop")
    Xtable <- category_combinations %>%
      left_join(Xtable_initial, by = cols) %>%
      mutate(n = case_when(
        is.na(n) ~ 0,
        !is.na(n) ~ n
      )) %>%
      as.data.frame() # sum of weights for treatments and untreated groups
    colnames(Xtable) <- c(colnames(Xmat), "n")
    SXY <- t(Y * weights) %*% Xmat # weighted sum of intercepts and treatments
    AD[["1"]] <- list(SXY = SXY, Xtable = Xtable)
  }
  AD
}

#' @keywords internal
computeWeights <- function(population, estimator = "ate") {
  if (estimator == "ate") {
    # 'Stabilized' ATE:
    return(ifelse(population$treatment == 1,
      mean(population$treatment == 1) / population$propensityScore,
      mean(population$treatment == 0) / (1 - population$propensityScore)
    ))
  } else {
    # 'Stabilized' ATT:
    return(ifelse(population$treatment == 1,
      mean(population$treatment == 1),
      mean(population$treatment == 0) * population$propensityScore / (1 - population$propensityScore)
    ))
  }
}

#' @keywords internal
# Function to fit GLM using aggregated data
#    allows for site-specific intercept by using heter_intercept = TRUE
oneshot_IPWGLM <- function(KSiteAD, formula, heter_intercept = FALSE) {
    K <- length(KSiteAD)
    if (heter_intercept == FALSE) {
        SXY <- KSiteAD[[1]]$SXY
        Xcat <- KSiteAD[[1]]$Xtable[, colnames(KSiteAD[[1]]$Xtable) != "n"]
        Xcat <- as.matrix(Xcat)
        counts <- KSiteAD[[1]]$Xtable$n
        if (K >= 2) {
            for (k in 2:K) {
                SXY <- SXY + KSiteAD[[k]]$SXY
                counts <- counts + KSiteAD[[k]]$Xtable$n
            }
        }
    } else {
        SXY.intercept <- KSiteAD[[1]]$SXY[1]
        SXY.cov <- KSiteAD[[1]]$SXY[-1]
        Xcat0 <- KSiteAD[[1]]$Xtable[, colnames(KSiteAD[[1]]$Xtable) != "n"]
        Xcat <- Xcat1 <- as.matrix(Xcat0[, -1])
        counts <- KSiteAD[[1]]$Xtable$n
        if (K >= 2) {
            for (k in 2:K) {
                SXY.intercept <- c(SXY.intercept, KSiteAD[[k]]$SXY[1])
                SXY.cov <- SXY.cov + KSiteAD[[k]]$SXY[-1]
                Xcat <- rbind(Xcat, Xcat1)
                counts <- c(counts, KSiteAD[[k]]$Xtable$n)
            }
        }
        SXY <- c(SXY.intercept, SXY.cov)
        SiteID <- rep(1:K, each = nrow(Xcat0))
        new.siteID <- sapply(1:K, function(i) ifelse(SiteID == i, 1, 0))
        colnames(new.siteID) <- paste0("Site", 1:K)
        Xcat <- cbind(new.siteID, Xcat)
        Xcat <- as.matrix(Xcat)
    }
    logLik_AD <- function(beta) {
        -(sum(SXY * beta) - sum(log(1 + exp(Xcat %*% c(beta))) * counts)) / sum(counts)
    }
    fit.AD <- optim(par = rep(0, ncol(Xcat)), logLik_AD, method = "BFGS")

    se <- sqrt(diag(solve(hessian(func = function(x) logLik_AD(x) * sum(counts), x = fit.AD$par))))
    # rownames(res) <- colnames(Xcat)
    confint <- fit.AD$par + c(-1, 1) * 1.96 * se
    # Store results in a list
    list(
      coefficients = fit.AD$par[2], 
      se = se[2],
      effect_size = exp(fit.AD$par),
      ci_lower = exp(confint[1]),
      ci_upper = exp(confint[2]),
      log_likelihood = -fit.AD$value * sum(counts),
      convergence = fit.AD$convergence
    ) 
}


################################################################################################################################
###################################################     Matching      ##########################################################
################################################################################################################################ 


#' @keywords internal
getAD_IPW_overlap <- function(SiteIPD, outcome_name, formula, link = "canonical", cut_off = NULL) {
  AD = list()
  Xmat <- grab_design_matrix(data = SiteIPD, rhs_formula = formula)
  Y <- SiteIPD[[outcome_name]]
  weights <- SiteIPD$weights
  Xmat.tbl <- data.frame(Xmat)
  category_combinations <- expand.grid(lapply(Xmat.tbl, unique),
    stringsAsFactors = FALSE
  ) %>% arrange_all()
  colnames(Xmat.tbl) <- colnames(category_combinations)

  if (link == "canonical") {
    Xmat.tbl <- as_tibble(Xmat.tbl)
    cols <- colnames(Xmat.tbl)
    Xmat.tbl.weights <- cbind(Xmat.tbl, weights = weights)
    Xtable_initial <- Xmat.tbl.weights %>%
      group_by_at(.vars = cols) %>%
      summarise(n = sum(weights), .groups = "drop")
    Xtable <- category_combinations %>%
      left_join(Xtable_initial, by = cols) %>%
      mutate(n = case_when(
        is.na(n) ~ 0,
        !is.na(n) ~ n
      )) %>%
      as.data.frame() # sum of weights for treatments and untreated groups
    colnames(Xtable) <- c(colnames(Xmat), "n")
    SXY <- t(Y * weights) %*% Xmat # weighted sum of intercepts and treatments
    AD[["1"]] <- list(SXY = SXY, Xtable = Xtable)
  }
  AD
}

#' @keywords internal
### overlapping weights 
computeWeights_overlap <- function(population, estimator = "ato") {
  if (estimator == "ato") {
    # 'Stabilized' ATE:
    return(ifelse(population$treatment == 1,
      mean(population$treatment == 1) * (1 - population$propensityScore),
      mean(population$treatment == 0) * (population$propensityScore)
    ))
  } else {
    # 'Stabilized' ATT:
    return(ifelse(population$treatment == 1,
      mean(population$treatment == 1),
      mean(population$treatment == 0) * population$propensityScore / (1 - population$propensityScore)
    ))
  }
}

################################################################################################################################
###########################     Other helper functions     ############################################
################################################################################################################################


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

#' @keywords internal
get_stratified_pop = function(mydata_test, nstrata){
  rowId = c(1:length(mydata_test$treatment))
  mydata_test <- cbind(rowId, mydata_test)
  stratifiedPop <- stratifyByPs(mydata_test, numberOfStrata = nstrata)
  return (stratifiedPop)
}

#' @keywords internal
stratifyByPs <- function(population, numberOfStrata = 5, stratificationColumns = c(), baseSelection = "all") {
  if (!("rowId" %in% colnames(population)))
    stop("Missing column rowId in population")
  if (!("treatment" %in% colnames(population)))
    stop("Missing column treatment in population")
  if (!("propensityScore" %in% colnames(population)))
    stop("Missing column propensityScore in population")
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

#' @keywords internal
get_SMD <-function(stratifiedPop, xvars){
  smd_before <- GetSMD(stratifiedPop[, xvars], stratifiedPop$treatment)
  weight_mat=Compute_weight(stratifiedPop)
  stratifiedPop=stratifiedPop%>%left_join(weight_mat,by="rowId")
  smd_after <- GetSMD(stratifiedPop[, xvars], stratifiedPop$treatment,stratifiedPop$weight)
  
  return(list(smd_before=smd_before,smd_after=smd_after))
}

#' @keywords internal
GetSMD <- function(data, treat, weights = NULL, std = TRUE){
  table <- col_w_smd(data, treat, weights, std)
  smd = t(data.frame(as.list(table)))
  smd = data.frame(row.names(smd), smd)
  rownames(smd) <- NULL
  colnames(smd) = c("Xvars", "SMD")
  
  return(smd)
}

#' @keywords internal
Compute_weight <- function(data) {
  stratumSize <- data %>%
    group_by(stratumId, treatment) %>%
    count() %>%
    ungroup()

  w <- stratumSize %>%
    mutate(weight = 1 / n) %>%
    inner_join(data, by = c("stratumId", "treatment"), multiple = "all") %>%
    dplyr::select(rowId, treatment, weight)

  wSum <- w %>%
    group_by(treatment) %>%
    summarize(wSum = sum(weight, na.rm = TRUE)) %>%
    ungroup()

  w_final <- w %>%
    inner_join(wSum, by = "treatment") %>%
    mutate(weight = weight / wSum) %>%
    dplyr::select(rowId, treatment, weight)

  return(w_final[, c(1, 3)])
}

#' Function to perform all data processing and pooled stratified analysis
#'
#' @param data The combined LATTE_ADRD data.
#' @param outcome_id The name of the primary outcome column.
#' @param outcome_time The name of the outcome time column (for Poisson).
#' @param sites A vector of site identifiers.
#' @return A list containing the results of the standard logistic and Poisson pooled analysis.
run_pooled_analysis <- function(data, outcome_id, outcome_time, sites) {
  
  n_sites <- length(sites)
  all_stratified_data <- NULL
  
  # Process each simulated site to calculate PS and strata
  for (site_id in 1:n_sites) {
    site_data <- data[data$site == site_id, ]
    
    # Identify covariates (xvars)
    xvars <- colnames(site_data)[!grepl("^outcome_", colnames(site_data)) & 
                                !colnames(site_data) %in% c("ID", "treatment", "index_date", "site", "group")]
    xvars <- xvars[colSums(site_data[xvars]) > 30]
    yvars <- colnames(site_data)[grepl("^outcome_", colnames(site_data))]
    
    # Prepare data for PS calculation
    mydata_ps <- site_data[, colnames(site_data) %in% c(xvars, "treatment", yvars, outcome_id, outcome_time)]
    
    # Calculate Propensity Scores (PS) using Lasso/Elastic Net
    form <- as.formula(paste("treatment ~ ", paste(xvars, collapse = "+")))
    Xmat <- grab_design_matrix(data = mydata_ps, rhs_formula = form)
    Y <- mydata_ps$treatment
    
    nfolds <- 10
    set.seed(42)
    foldid <- sample(rep(seq(nfolds), length.out = length(Y)))
    
    Fit_ps_cv <- tryCatch({
      cv.glmnet(Xmat, Y, alpha = 1, family = "binomial", nfolds = nfolds, foldid = foldid)
    }, error = function(e) {
      message("CV failed for site ", site_id, ". Using fixed lambda.")
      list(lambda.min = 0.01)
    })
    Fit_ps <- glmnet(Xmat, Y, alpha = 1, family = "binomial", lambda = Fit_ps_cv$lambda.min)
    propensityScore <- predict(Fit_ps, (Xmat), type = "response")
    mydata_ps$propensityScore <- propensityScore
    
    # Create stratified population
    best_strata <- optimize_strata(mydata_ps, xvars)
    stratifiedPop <- get_stratified_pop(mydata_test = mydata_ps, nstrata = best_strata$n_strata)
    
    # Add identifiers
    stratifiedPop$site <- site_id
    stratifiedPop$global_stratumId <- paste0(site_id, "_", stratifiedPop$stratumId)
    
    # Collect all stratified data
    if (is.null(all_stratified_data)) {
      all_stratified_data <- stratifiedPop
    } else {
      common_cols <- intersect(colnames(all_stratified_data), colnames(stratifiedPop))
      all_stratified_data <- rbind(
        all_stratified_data[, common_cols],
        stratifiedPop[, common_cols]
      )
    }
  }
  
  # --- Standard Pooled Analysis ---
  
  # Logistic Regression
  outcome_formula <- as.formula(paste0(outcome_id, "~ treatment + factor(global_stratumId)"))
  fit_log <- glm(outcome_formula, data = all_stratified_data, family = binomial(link = "logit"))
  
  log_res <- summary(fit_log)$coefficients["treatment", ]
  normal_est <- log_res["Estimate"]
  normal_se <- log_res["Std. Error"]
  
  
  return(list(
    logistic = list(est = normal_est, se = normal_se)
  ))
}

#' @keywords internal
## Function to neatly print the results
print_results <- function(title, est, se) {
  effect_size <- exp(est)
  ci_lower <- exp(est - 1.96 * se)
  ci_upper <- exp(est + 1.96 * se)
  
  cat(paste("\n", title, ":\n", sep = ""))
  cat("Coefficient:", round(est, 2), "\n")
  cat("Standard Error:", round(se, 2), "\n")
  cat("Odds Ratio/Rate Ratio:", round(effect_size, 2), "\n")
  cat("95% CI:", sprintf("[%.2f, %.2f]", ci_lower, ci_upper), "\n")
}
