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
  xvars = control$variables
  # Create PS model formula
  ps_formula <- as.formula(paste("treatment ~", paste(xvars, collapse = "+")))
  outcomes = c(control$outcome, control$nco_outcomes)
  yvars <- outcomes
  if (control$outcome_model == "poisson") { 
    outcomes_time = control$nco_outcome_times
    yvars <- c(outcomes, outcomes_time)
  }
  ADdatas = list()

  # Prepare data for PS model
  

  mydata_ps <- ipdata[, colnames(ipdata) %in% c(xvars, "treatment", yvars)]

  Xmat <- grab_design_matrix(data = mydata_ps, rhs_formula = ps_formula)
  Y <- mydata_ps$treatment
  
  # Fit PS model using cv.glmnet
  nfolds <- 10
  set.seed(42)
  foldid <- sample(rep(seq(nfolds), length.out = length(Y)))
  
  ps_fit_cv <- cv.glmnet(Xmat, Y, alpha = 1, family = "binomial", 
                        nfolds = nfolds, foldid = foldid)
  ps_fit <- glmnet(Xmat, Y, alpha = 1, family = "binomial", 
                  lambda = ps_fit_cv$lambda.min)
    
  # Calculate propensity scores
  propensityScore <- predict(ps_fit, Xmat, type = "response")
  mydata_ps$propensityScore <- propensityScore
  for (i in 1:length(outcomes)) {
    outcome = outcomes[i]
    if (control$outcome_model == "poisson") {
      outcome_time = outcomes_time[i]
    }
    
    outcome_formula <- as.formula(paste0( outcome, " ~ treatment"))
    if (control$balancing_method == "stratification") {
      # Find optimal number of strata
      best_strata <- optimize_strata(mydata_ps, xvars)
      stratifiedPop <- get_stratified_pop(mydata_ps, nstrata = best_strata$n_strata)
      if (control$outcome_model == "logistic") {
        ADdata <- create_2x2_tables(
          stratifiedPop$treatment,
          stratifiedPop[[outcome]],
          stratifiedPop$stratumId
        )
      } else if (control$outcome_model == "poisson") {
        ADdata <- create_2x2_tables_poisson(
          stratifiedPop,
          outcome_id = outcome,
          outcome_time = outcome_time
        )
      }
    } else if (control$balancing_method == "matching") {
      IPTW_ato <- computeWeights_overlap(mydata_ps)
      mydata_ps$weights <- IPTW_ato
      ADdata <- getAD_IPW_overlap(mydata_ps, paste0(outcome_id), outcome_formula)
    } else if (control$balancing_method == "IPTW") {
      result <- trimByW(mydata_ps$propensityScore, mydata_ps$treatment, trimFraction = 0.05)
      mydata_ps <- mydata_ps[result, ]
      IPTW <- computeWeights(mydata_ps)
      mydata_ps$weights <- IPTW
      ADdata <- getAD_IPW(mydata_ps, paste0(outcome_id), outcome_formula)
    } else {
      stop("Unknown PS method specified in config. Must be one of 'stratification', 'matching', or 'IPTW'")
    }
    ADdatas[[outcome]] <- ADdata
  }
  # Return initialization results
  return(list(
    prepared_data = ADdatas
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
LATTE.estimate <- function(init_data, control, config) { 
  K <- length(control$sites)
  results_all = list()
  outcomes = c(control$outcome, control$nco_outcomes)
  results_table <- data.frame(
    outcome = character(),
    est = numeric(),   # log effect (log OR or log RR)
    se  = numeric(),
    is_nco = logical(),
    stringsAsFactors = FALSE
  )
  for (outcome in outcomes){
    ADdata <- list()
    results <- list()
    for (site_i in control$sites) {
      init_i <- pdaGet(paste0(site_i, "_initialize"), config)
      outcome_data = init_i$prepared_data[[outcome]]
      if (length(outcome_data) != 0) {
        for (j in 1:length(outcome_data)) {
          ADdata[[length(ADdata) + 1]] <- outcome_data[[j]]
        }
      }
    }
    if (length(ADdata) == 0) {
      next
    }
    outcome_formula <- as.formula(paste0(outcome, " ~ treatment"))
    # Fit conditional logistic regression
    if (control$balancing_method == "stratification") {
      if (control$outcome_model == "logistic") {
        results <- optimize_conditional_logistic_2x2(ADdata)
      } else if (control$outcome_model == "poisson") {
        results <- optimize_conditional_poisson_2x2(ADdata)
      }
    } else if (control$balancing_method == "matching") {
      results <- oneshot_IPWGLM(ADdata, outcome_formula)
    } else if (control$balancing_method == "IPTW") {
      results <- oneshot_IPWGLM(ADdata, outcome_formula)
    }
    results_all[[outcome]] <- list(
      coefficients = results$coefficients,
      se = results$se,
      effect_size = results$effect_size,
      ci_lower = results$ci_lower,
      ci_upper = results$ci_upper,
      convergence = results$convergence
    )
    # Append to tidy table for calibration pool
    # If your helpers expose the treatment coefficient differently,
    # replace results$effect_size and results$se below accordingly.
    results_table <- rbind(
      results_table,
      data.frame(
        outcome = outcome,
        est = as.numeric(results$coefficients),
        se  = as.numeric(results$se),
        is_nco = outcome %in% control$nco_outcomes,
        stringsAsFactors = FALSE
      )
    )
  }
  # }
  # Prepare negative control set
  neg_df <- subset(
    results_table,
    is_nco & !is.na(est) & !is.na(se) & is.finite(est) & is.finite(se) & abs(est) <= 10
  )
  
  # Fit empirical null (for p-value calibration) and systematic error model (for CI calibration)
  fitnull <- fitNull(logR = neg_df$est, seLogRr = neg_df$se)
  se_model <- fitSystematicErrorModel(neg_df$est,
                                      neg_df$se,
                                      rep(0, nrow(neg_df)))  # NCO truth = 0 on log scale
  # Apply calibration to *non*-NCO outcomes
  target_df <- subset(results_table, !is_nco & !is.na(est) & !is.na(se) & is.finite(est) & is.finite(se))
  if (nrow(target_df) > 0) {
    # Calibrated CIs
    cal_ci <- calibrateConfidenceInterval(
      logRr = target_df$est,
      seLogRr = target_df$se,
      model = se_model,
      ciWidth = 0.95
    )
    colnames(cal_ci) <- c("cal_est", "cal_ll", "cal_ul", "cal_se")
    # Calibrated p-values
    cal_p <- calibrateP(null = fitnull, logR = target_df$est, seLogRr = target_df$se)
    # Merge calibrated stats back into the per-outcome list
    for (i in seq_len(nrow(target_df))) {
      outname <- target_df$outcome[i]
      results_all[[outname]]$calibrated <- list(
        est = as.numeric(cal_ci$cal_est[i]),
        se  = as.numeric(cal_ci$cal_se[i]),
        effect_size = exp(as.numeric(cal_ci$cal_est[i])),
        ci_lower = exp(as.numeric(cal_ci$cal_ll[i])),
        ci_upper = exp(as.numeric(cal_ci$cal_ul[i])),
        p_value = as.numeric(cal_p[i]),
        ci_width = 0.95
      )
    }
  }
  # Also annotate NCO entries with their (uncalibrated) traditional p and calibrated p for QC
  nco_trad_p <- computeTraditionalP(logR = neg_df$est, seLogRr = neg_df$se)
  nco_cal_p  <- calibrateP(null = fitnull, logR = neg_df$est, seLogRr = neg_df$se)
  for (i in seq_len(nrow(neg_df))) {
    outname <- neg_df$outcome[i]
    results_all[[outname]]$nco_qc <- list(
      traditional_p = as.numeric(nco_trad_p[i]),
      calibrated_p  = as.numeric(nco_cal_p[i])
    )
  }
  
  return(list(
    by_outcome = results_all,
    calibration = list(
      used_nco_count = nrow(neg_df)
    ),
    raw_table = results_table
  ))
}
 
