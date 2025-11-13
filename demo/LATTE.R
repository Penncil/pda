# Demo script for LATTE
# Load required packages
# Assuming 'pda' package is loaded or source files are available for pda() and getCloudConfig()
# source("R/pda.R")
# source("R/LATTE.R")
# source("R/latte-misc.R")
## In the toy example below we aim to analyze the treatment effects of acetaminophen on ADRD using logistic regression, and propensity score stratification,
## data: latte_synthetic_data.rda, we randomly assign to 3 sites: 'site1', 'site2', 'site3'
## we demonstrate using PDA LATTE can obtain a surrogate estimator that is close to the pooled estimate.
## We run the example in local directory. In actual collaboration, account/password for pda server
## will be assigned to the sites at the server https://pda.one.
## Each site can access via web browser to check the communication of the summary stats.


# Define variables and load data
set.seed(42)
outcome_id = "outcome_ADRD_value"
outcome_time = "outcome_ADRD_time"
sites <- c("site1", "site2", "site3")
data("LATTE_ADRD", package = "pda")

# Separate data for LATTE analysis
site_data <- list(
  site1 = LATTE_ADRD[LATTE_ADRD$site == 1, ],
  site2 = LATTE_ADRD[LATTE_ADRD$site == 2, ],
  site3 = LATTE_ADRD[LATTE_ADRD$site == 3, ]
)

# Identify covariates from site1 data for control list
site1_data_subset <- LATTE_ADRD[LATTE_ADRD$site == 3, ]
xvars <- colnames(site1_data_subset)[!grepl("^outcome_", colnames(site1_data_subset)) & 
                            !colnames(site1_data_subset) %in% c("ID", "treatment", "index_date", "site", "group")]
xvars <- xvars[colSums(site1_data_subset[xvars]) > 30]

# Non-Concurrent Outcomes (NCOs) setup for LATTE calibration
nco_outcomes = c(
  "acute_conjunctivitis", "acute_tonsillitis", "adhesive_capsulitis_of_shoulder", "allergic_rhinitis", 
  "blepharitis", "carpal_tunnel_syndrome", "chalazion", "contact_dermatitis", "dental_caries", 
  "deviated_nasal_septum", "foreign_body_in_ear", "gout", "hemorrhoids", "impacted_cerumen", 
  "influenza", "ingrowing_nail", "low_back_pain", "menieres_disease", "osteoarthritis_of_knee", 
  "osteoporosis", "foot_drop", "hearing_problem", "intra_abdominal_and_pelvic_swelling_mass_and_lump", 
  "irritability_and_anger", "wristdrop"
)
nco_outcomes <- paste0("outcome_", nco_outcomes, "_value")
outcome_times = c(outcome_time, paste0("outcome_", c("acute_conjunctivitis", "acute_tonsillitis", 
                                                    "adhesive_capsulitis_of_shoulder", "allergic_rhinitis", 
                                                    "blepharitis", "carpal_tunnel_syndrome", "chalazion", 
                                                    "contact_dermatitis", "dental_caries", 
                                                    "deviated_nasal_septum", "foreign_body_in_ear", 
                                                    "gout", "hemorrhoids", "impacted_cerumen", 
                                                    "influenza", "ingrowing_nail", "low_back_pain", 
                                                    "menieres_disease", "osteoarthritis_of_knee", 
                                                    "osteoporosis", "foot_drop", "hearing_problem", 
                                                    "intra_abdominal_and_pelvic_swelling_mass_and_lump", 
                                                    "irritability_and_anger", "wristdrop"), "_time"))


# --- Setup: Working directory ---
if (!dir.exists("pda_latte_results")) {
  dir.create("pda_latte_results")
}
original_wd <- getwd()
setwd("pda_latte_results")

# 1. Run Pooled Analysis (Traditional Benchmark)
cat("## Running Traditional Pooled Analysis...\n")
pooled_results <- run_pooled_analysis(LATTE_ADRD, outcome_id, outcome_time, sites)

# 2. Setup LATTE Control Parameters
control <- list(
  project_name = "LATTE Demo Study",
  step = "initialize", # Will be set by run_latte_analysis
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
  balancing_method = "stratification", # stratification
  outcome_model = "logistic",          # logistic
  nco_outcome_times = outcome_times    # needed if outcome_model is poisson
)

# 3. Run LATTE Analysis (Distributed Simulation)
cat("\n## Running Distributed LATTE Analysis (Initialize & Estimate)...\n")
dir.create("pda_latte_results")
setwd("pda_latte_results")
# Step 1: Initialize at each site (each site computes their 2x2 tables)
menu <- function(choices, title = NULL) 1

pda(site_id = 'site1', control = control, ipdata = site_data$site1, dir = getwd())

menu <- function(choices, title = NULL) 1

pda(site_id = 'site2', ipdata = site_data$site2, dir = getwd())
menu <- function(choices, title = NULL) 1

pda(site_id = 'site3', ipdata = site_data$site3, dir = getwd())
menu <- function(choices, title = NULL) 1

pda(site_id = "site1", ipdata = site_data$site1, dir = getwd())

#############################  STEP 2: estimate  ###############################
menu <- function(choices, title = NULL) 1

pda(site_id = "site1", ipdata = site_data$site1, dir = getwd())


# Get the results using pdaGet
config <- getCloudConfig(site_id = 'site1', dir = getwd())
latte_results <- pdaGet('site1_estimate', config)

# --- Results Comparison ---
cat("\n\n######################################################################")
cat("\n#                 Comparison of Results: Acetaminophen vs ADRD         #")
cat("\n######################################################################\n")
# Pooled Logistic
print_results("Pooled Stratified Logistic Regression", 
              pooled_results$logistic$est, pooled_results$logistic$se)

cat("\n---\n")

# LATTE Uncalibrated
print_results("LATTE Distributed Analysis (Uncalibrated)",
              latte_results$by_outcome[[outcome_id]]$coefficients,
              latte_results$by_outcome[[outcome_id]]$se)

cat("\n---\n")

# LATTE Calibrated
print_results("LATTE Distributed Analysis (NCO Calibrated)",
              latte_results$by_outcome[[outcome_id]]$calibrated$est,
              latte_results$by_outcome[[outcome_id]]$calibrated$se)

# Clean up
setwd(original_wd)
cat(paste0("\n\nResults files are located in the '", original_wd, "/pda_latte_results' directory.\n"))
