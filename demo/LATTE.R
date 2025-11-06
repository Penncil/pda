# Demo script for LATTE 
# Load required packages
require(survival)
require(data.table)
require(pda)
require(glmnet)
require(dplyr)
require(Matrix)
require(tibble)
library(cobalt)
library(geex)
source("/Users/luli/pda2/pda/R/pda.R")
source("/Users/luli/pda2/pda/R/LATTE.R")
source("/Users/luli/pda2/pda/LATTE_codes/latte_codes.R")

# Create sample data
set.seed(42)
outcome_id = "ADRD"
n <- 8000  # total sample size 
sites <- c("site1", "site2", "site3")
n_sites = length(sites)
cohort = read.csv("/Users/luli/pda2/pda/161all_data_sim.csv")
cohort <- cohort %>% select(-X)
cohort <- cohort %>%
    mutate(across(ends_with("_value"), ~ replace(., . == -1, 0)))

cohort <- cohort %>%
    mutate(site = sample(1:n_sites, n(), replace = TRUE))

xvars <- colnames(cohort)[!grepl("^outcome_", colnames(cohort)) & 
                               !colnames(cohort) %in% c("ID", "treatment", "index_date", "site", "group")]
# Arrays to store data for standard logistic regression
all_stratified_data <- NULL
KSiteAD_uf <- list()

# Process each simulated site
for (site_id in 1:n_sites) {
  # Extract data for this site
  site_data <- cohort[cohort$site == site_id, ]
  
  # Process covariates
  xvars <- colnames(site_data)[!grepl("^outcome_", colnames(site_data)) & 
                              !colnames(site_data) %in% c("ID", "treatment", "index_date", "site", "group")]
  xvars <- xvars[colSums(site_data[xvars]) > 30]
  yvars <- colnames(site_data)[grepl("^outcome_", colnames(site_data))]
  
  # Calculate propensity scores for this site
  form <- as.formula(paste("treatment ~ ", paste(xvars, collapse = "+")))
  mydata_ps <- site_data[, colnames(site_data) %in% c(xvars, "treatment", yvars)]
  
  # Add the outcome variable if not already included
  outcome_var <- paste0("outcome_", outcome_id, "_value")
  if (!(outcome_var %in% colnames(mydata_ps)) && outcome_var %in% colnames(site_data)) {
    mydata_ps[[outcome_var]] <- site_data[[outcome_var]]
  }
  
  Xmat <- grab_design_matrix(data = mydata_ps, rhs_formula = form)
  Y <- mydata_ps$treatment
  
  nfolds <- 10
  foldid <- sample(rep(seq(nfolds), length.out = length(Y)))
  set.seed(42)
  Fit_ps_cv <- tryCatch({
    cv.glmnet(Xmat, Y, alpha = 1, family = "binomial", nfolds = nfolds, foldid = foldid)
  }, error = function(e) {
    # If cross-validation fails, use a small fixed lambda
    message("CV failed for site ", site_id, ". Using fixed lambda.")
    list(lambda.min = 0.01)
  })
  
  Fit_ps <- glmnet(Xmat, Y, alpha = 1, family = "binomial", lambda = Fit_ps_cv$lambda.min)
  propensityScore <- predict(Fit_ps, (Xmat), type = "response")
  mydata_ps$propensityScore <- propensityScore
  
  # Create stratified population for this site
  best_strata <- optimize_strata(mydata_ps, xvars)
  stratifiedPop <- get_stratified_pop(mydata_test = mydata_ps, nstrata = best_strata$n_strata)
  
  # Add site identifier to the stratified data
  stratifiedPop$site <- site_id
  # Create unique stratum IDs across all sites
  stratifiedPop$global_stratumId <- paste0(site_id, "_", stratifiedPop$stratumId)
  
  # Collect all stratified data for the standard logistic regression
  if (is.null(all_stratified_data)) {
    all_stratified_data <- stratifiedPop
  } else {
    # Ensure columns match before rbind
    common_cols <- intersect(colnames(all_stratified_data), colnames(stratifiedPop))
    all_stratified_data <- rbind(
      all_stratified_data[, common_cols],
      stratifiedPop[, common_cols]
    )
  }
  # Create 2x2 tables for this site for conditional logistic regression
    for (strat in unique(stratifiedPop$stratumId)) {
      strat_data <- stratifiedPop[stratifiedPop$stratumId == strat, ]
      
      # Only include strata with variation in treatment and outcome
      if (length(unique(strat_data$treatment)) > 1 && 
          length(unique(strat_data[[paste0("outcome_", outcome_id, "_value")]])) > 1) {
        
        # Create the 2x2 contingency table
        table_2x2 <- matrix(0, nrow = 2, ncol = 2)
        
        # Exposed cases (treatment=1, outcome=1)
        table_2x2[1,1] <- sum(strat_data$treatment == 1 & 
                               strat_data[[paste0("outcome_", outcome_id, "_value")]] == 1)
        
        # Unexposed cases (treatment=0, outcome=1)
        table_2x2[1,2] <- sum(strat_data$treatment == 0 & 
                               strat_data[[paste0("outcome_", outcome_id, "_value")]] == 1)
        
        # Exposed controls (treatment=1, outcome=0)
        table_2x2[2,1] <- sum(strat_data$treatment == 1 & 
                               strat_data[[paste0("outcome_", outcome_id, "_value")]] == 0)
        
        # Unexposed controls (treatment=0, outcome=0)
        table_2x2[2,2] <- sum(strat_data$treatment == 0 & 
                               strat_data[[paste0("outcome_", outcome_id, "_value")]] == 0)
        
        # Only add tables that are informative
        if (sum(table_2x2) > 0 && min(rowSums(table_2x2)) > 0 && min(colSums(table_2x2)) > 0) {
          KSiteAD_uf[[length(KSiteAD_uf) + 1]] <- table_2x2
        }
      }
    }
}
# Method 1: Standard stratified logistic regression with site and stratum as fixed effects
outcome_formula <- as.formula(paste0("outcome_", outcome_id, "_value ~ treatment + factor(global_stratumId)"))
# Check if we have the necessary columns
required_cols <- c(paste0("outcome_", outcome_id, "_value"), "treatment", "global_stratumId")
missing_cols <- setdiff(required_cols, colnames(all_stratified_data))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse=", "))
}
unique(all_stratified_data$global_stratumId)
# Run the standard logistic regression
fit <- tryCatch({
  glm(outcome_formula, data = all_stratified_data, family = binomial(link = "logit"))
}, error = function(e) {
  message("Error in standard logistic regression: ", e$message)
  # Return a dummy model with NA coefficients
  structure(list(
    coefficients = c("(Intercept)" = NA, "treatment" = NA),
    rank = 2,
    family = binomial(link = "logit"),
    linear.predictors = rep(NA, nrow(all_stratified_data)),
    fitted.values = rep(NA, nrow(all_stratified_data)),
    residuals = rep(NA, nrow(all_stratified_data)),
    df.residual = nrow(all_stratified_data) - 2,
    converged = FALSE
  ), class = "glm")
})

# Extract standard logistic regression results
if (!is.na(fit$coefficients["treatment"])) {
  normal_est <- summary(fit)$coefficients["treatment", "Estimate"]
  normal_se <- summary(fit)$coefficients["treatment", "Std. Error"]
  normal_pval <- summary(fit)$coefficients["treatment", "Pr(>|z|)"]
} 

ResAD <- optimize_conditional_logistic_2x2(KSiteAD_uf)
cond_est <- ResAD$coefficients
cond_se <- ResAD$se
cond_pval <- 2 * (1 - pnorm(abs(cond_est / cond_se)))


site_data <- list(
  site1 = cohort[cohort$site == 1, ],
  site2 = cohort[cohort$site == 2, ],
  site3 = cohort[cohort$site == 3, ]
) 

  
# Initialize LATTE analysis
control <- list(
  project_name = "LATTE Demo Study",
  step = "initialize",
  sites = sites,
  model = "LATTE",
  family = "binomial",
  outcome = "outcome_ADRD_value",
  variables = xvars,
  lead_site = "site1",
  min_count = 5,
  nfolds = 5,
  min_strata = 2,
  max_strata = 6,
  start_beta = 0.1,
  ### option for choosing method
  balancing_method = "stratification" # overlapping, IPTW
)
set.seed(42)
setwd("/Users/luli/pda2/pda/test")
# Step 1: Initialize at each site (each site computes their 2x2 tables)
pda(site_id = 'site1', control = control, ipdata = site_data$site1, dir = getwd())

pda(site_id = 'site2', ipdata = site_data$site2, dir = getwd())
pda(site_id = 'site3', ipdata = site_data$site3, dir = getwd())

# auto change to estimate 
pda(site_id = "site1", ipdata = site_data$site1, dir = getwd())


# Get the results using pdaGet
config <- getCloudConfig(site_id = 'site1', dir = getwd())
cat("Working directory:", getwd(), "\n")
cat("Config contents:", str(config), "\n")
config$dir = "/Users/luli/pda2/pda/test"
config$dir <- gsub("\\\\", "/", config$dir)

latte_results <- pdaGet('site1_estimate', config)

# Print comparison of results
cat("\nComparison of Results:\n")
cat("\nTrue treatment effect: 0.3\n")





cat("\nPooled Analysis (Traditional):\n")
print(summary(fit.pool)$coefficients["treatment", ])
exp(summary(fit.pool)$coefficients["treatment", ])

cat("\nLATTE Analysis:\n")
cat("Coefficient:", latte_results$coefficients, "\n")
cat("Standard Error:", latte_results$se, "\n")
cat("Odds Ratio:", latte_results$odds_ratio, "\n")
cat("95% CI:", sprintf("[%f, %f]", latte_results$ci_lower, latte_results$ci_upper), "\n")

# Optional: Print convergence information
cat("\nConvergence Information:\n")
cat("Convergence status:", latte_results$convergence, "\n")
cat("Message:", latte_results$message, "\n")

