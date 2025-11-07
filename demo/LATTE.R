# Install devtools if not installed
# install.packages("devtools")

# Install vscDebugger from GitHub
# devtools::install_github("ManuelHentschel/vscDebugger")
library(cli)
# install.packages(c("pillar","rlang","lifecycle"))
# install.packages("remotes", repos = "https://cloud.r-project.org")
# remotes::install_version("cli", version = "3.6.5", repos = "https://cloud.r-project.org")
# packageVersion("cli")

# Demo script for LATTE 
# Load required packages
require(vscDebugger)
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
  time_var <- paste0("outcome_", outcome_id, "_time")
  if (!(outcome_var %in% colnames(mydata_ps)) && outcome_var %in% colnames(site_data)) {
    mydata_ps[[outcome_var]] <- site_data[[time_var]]
  }
  
  Xmat <- grab_design_matrix(data = mydata_ps, rhs_formula = form)
  Y <- mydata_ps$treatment
  
  nfolds <- 10
  set.seed(42)
  foldid <- sample(rep(seq(nfolds), length.out = length(Y)))
  
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
}
# Method 1: Standard stratified logistic regression with site and stratum as fixed effects
outcome_formula <- as.formula(paste0("outcome_", outcome_id, "_value ~ treatment + factor(global_stratumId)"))
# Check if we have the necessary columns
required_cols <- c(paste0("outcome_", outcome_id, "_value"), "treatment", "global_stratumId")
missing_cols <- setdiff(required_cols, colnames(all_stratified_data))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Run the standard logistic regression
fit <- glm(outcome_formula, data = all_stratified_data, family = binomial(link = "logit"))

# Method 2: Standard stratified Poisson regression with site and stratum as fixed effects
outcome_formula_poisson <- as.formula(paste0("outcome_", outcome_id, "_value ~ treatment + factor(global_stratumId) + offset(offset_term)"))
offset_term <- log(all_stratified_data[[time_var]])
fit_poisson <- glm(outcome_formula_poisson, data = all_stratified_data, family = poisson(link = "log"))

# Extract standard logistic regression results
if (!is.na(fit$coefficients["treatment"])) {
  normal_est <- summary(fit)$coefficients["treatment", "Estimate"]
  normal_se <- summary(fit)$coefficients["treatment", "Std. Error"]
  normal_pval <- summary(fit)$coefficients["treatment", "Pr(>|z|)"]
}
# Extract standard Poisson regression results
if (!is.na(fit_poisson$coefficients["treatment"])) {
  normal_est_poisson <- summary(fit_poisson)$coefficients["treatment", "Estimate"]
  normal_se_poisson <- summary(fit_poisson)$coefficients["treatment", "Std. Error"]
  normal_pval_poisson <- summary(fit_poisson)$coefficients["treatment", "Pr(>|z|)"]
}
 
# Prepare data for LATTE analysis
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
  balancing_method = "stratification", # overlapping, IPTW
  ### option for outcome model 
  outcome_model = "poisson", # logistic, poisson
  outcome_time = "outcome_ADRD_time" # only needed if outcome_model is poisson
)

setwd("/Users/luli/pda2/pda/test")
# Step 1: Initialize at each site (each site computes their 2x2 tables)
menu <- function(choices, title = NULL) 1

pda(site_id = 'site1', control = control, ipdata = site_data$site1, dir = getwd())

menu <- function(choices, title = NULL) 1

pda(site_id = 'site2', ipdata = site_data$site2, dir = getwd())
menu <- function(choices, title = NULL) 1

pda(site_id = 'site3', ipdata = site_data$site3, dir = getwd())
menu <- function(choices, title = NULL) 1

pda(site_id = "site1", ipdata = site_data$site1, dir = getwd())
# auto change to estimate
menu <- function(choices, title = NULL) 1

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

cat("\nPooled Analysis (Traditional):\n")
cat("Coefficient:", round(normal_est, 2), "\n")
cat("Standard Error:", round(normal_se, 2), "\n")
cat("Odds Ratio:", round(exp(normal_est), 2), "\n")
cat("95% CI:", sprintf("[%.2f, %.2f]", exp(normal_est - 1.96 * normal_se), exp(normal_est + 1.96 * normal_se)), "\n")



cat("\nLATTE Analysis:\n")
cat("Coefficient:", round(latte_results$coefficients, 2), "\n")
cat("Standard Error:", round(latte_results$se, 2), "\n")
cat("Odds Ratio:", round(latte_results$effect_size, 2), "\n")
cat("95% CI:", sprintf("[%.2f, %.2f]", latte_results$ci_lower, latte_results$ci_upper), "\n")

# # Optional: Print convergence information
# cat("\nConvergence Information:\n")
# cat("Convergence status:", latte_results$convergence, "\n")
# cat("Message:", latte_results$message, "\n")
 
