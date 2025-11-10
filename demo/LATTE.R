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
library(numDeriv)
source("R/pda.R")
source("R/LATTE.R")
source("R/LATTE_helper.R")
library(EmpiricalCalibration)

## In the toy example below we aim to analyze the treatment effects of acetaminophen on ADRD using logistic regression, and propensity score stratification,
## data: latte_synthetic_data.rda, we randomly assign to 3 sites: 'site1', 'site2', 'site3'
## we demonstrate using PDA LATTE can obtain a surrogate estimator that is close to the pooled estimate.
## We run the example in local directory. In actual collaboration, account/password for pda server
## will be assigned to the sites at the server https://pda.one.
## Each site can access via web browser to check the communication of the summary stats.

# define variables and load data
set.seed(42)

outcome_id = "outcome_ADRD_value"
nco_outcomes = c(
  "acute_conjunctivitis", "acute_tonsillitis", "adhesive_capsulitis_of_shoulder", "allergic_rhinitis", "blepharitis",
  "carpal_tunnel_syndrome", "chalazion", "contact_dermatitis", "dental_caries", "deviated_nasal_septum", "foreign_body_in_ear",
  "gout", "hemorrhoids", "impacted_cerumen", "influenza", "ingrowing_nail", "low_back_pain", "menieres_disease", "osteoarthritis_of_knee",
  "osteoporosis", "foot_drop", "hearing_problem", "intra_abdominal_and_pelvic_swelling_mass_and_lump", "irritability_and_anger",
  "wristdrop"
)
nco_outcomes_time = paste0("outcome_", nco_outcomes, "_time")

nco_outcomes <- paste0("outcome_", nco_outcomes, "_value")
outcome_time = "outcome_ADRD_time"
outcome_times = c(outcome_time, nco_outcomes_time)
sites <- c("site1", "site2", "site3")
n_sites = length(sites)

load("/Users/luli/pda2/pda/data/latte_synthetic_data.rda")
 

######################################################################
############ Get Pooled results      #################################
######################################################################

all_stratified_data <- NULL
KSiteAD_uf <- list()

# # Process each simulated site
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
  outcome_var <- outcome_id
  if (!(outcome_var %in% colnames(mydata_ps)) && outcome_var %in% colnames(site_data)) {
    mydata_ps[[outcome_var]] <- site_data[[outcome_var]]
  }
  time_var <- outcome_time
  if (!(time_var %in% colnames(mydata_ps)) && time_var %in% colnames(site_data)) {
    mydata_ps[[time_var]] <- site_data[[time_var]]
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
outcome_formula <- as.formula(paste0(outcome_id, "~ treatment + factor(global_stratumId)"))
# Check if we have the necessary columns
required_cols <- c(paste0( outcome_id ), "treatment", "global_stratumId")
missing_cols <- setdiff(required_cols, colnames(all_stratified_data))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Run the standard logistic regression
fit <- glm(outcome_formula, data = all_stratified_data, family = binomial(link = "logit"))

# Method 2: Standard stratified Poisson regression with site and stratum as fixed effects
outcome_formula_poisson <- as.formula(paste0(outcome_id, " ~ treatment + factor(global_stratumId) + offset(offset_term)"))
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

######################################################################
############ Run LATTE results      ##################################
######################################################################



# Initialize LATTE analysis
control <- list(
  project_name = "LATTE Demo Study",
  step = "initialize",
  sites = sites,
  model = "LATTE",
  family = "binomial",
  outcome = outcome_id,
  nco_outcomes = nco_outcomes,
  variables = xvars,
  lead_site = "site1",
  min_count = 5,
  nfolds = 5,
  min_strata = 2,
  max_strata = 6,
  start_beta = 0.1,
  ### option for choosing method
  balancing_method = "stratification", # matching, IPTW, stratification
  ### option for outcome model 
  outcome_model = "logistic", # logistic, poisson
  nco_outcome_times = outcome_times# only needed if outcome_model is poisson
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
cat("Coefficient:", round(latte_results$by_outcome[[outcome_id]]$coefficients, 2), "\n")
cat("Standard Error:", round(latte_results$by_outcome[[outcome_id]]$se, 2), "\n")
cat("Odds Ratio:", round(latte_results$by_outcome[[outcome_id]]$effect_size, 2), "\n")
cat("95% CI:", sprintf("[%.2f, %.2f]", latte_results$by_outcome[[outcome_id]]$ci_lower,latte_results$by_outcome[[outcome_id]]$ci_upper), "\n")
 
 
######################################################################
############ NCO calibrated results ##################################
######################################################################



cat("\nLATTE calibrated Analysis:\n")
cat("Coefficient:", round(latte_results$by_outcome[[outcome_id]]$calibrated$est, 2), "\n")
cat("Standard Error:", round(latte_results$by_outcome[[outcome_id]]$calibrated$se, 2), "\n")
cat("Odds Ratio:", round(latte_results$by_outcome[[outcome_id]]$calibrated$effect_size, 2), "\n")
cat("95% CI:", sprintf("[%.2f, %.2f]", latte_results$by_outcome[[outcome_id]]$calibrated$ci_lower, latte_results$by_outcome[[outcome_id]]$calibrated$ci_upper), "\n")
 
