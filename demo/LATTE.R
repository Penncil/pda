# Demo script for LATTE
# Load required packages
# Assuming 'pda' package is loaded or source files are available for pda() and getCloudConfig()
# source("R/pda.R")
# source("R/LATTE.R")
# source("R/latte-misc.R")
## In the toy example below we aim to analyze the treatment effects of acetaminophen on ADRD using logistic regression, and propensity score stratification,
## data: latte_synthetic_data.rda, we randomly assign to 3 sites: 'site1', 'site2', 'site3'
## we demonstrate using PDA LATTE can obtain a surrogate estimator that is close to the pooled estimate.
## We run the example in local directory.
## In actual collaboration, the data communication can be done via the PDA_OTA platform https://pda-ota.pdamethods.org/
## Each site can access via web browser to transfer aggregate data and check the progress of the project.
# Define variables and load data
set.seed(42)
outcome_id = "status"
outcome_time = "time"
treatment_var = "Trt"
xvars = c("Age", "Sex", "RE", "Mutation")

sites <- c("site1", "site2", "site3", "site4", "site5", "site6", "site7", "site8", "site9", "site10")

# data = read.csv("/Users/luli/Documents/developer/pda1210/pda/JJ_pda_simu_data_20251123.csv")
data = read.csv("../../JJ_pda_simu_data_20251123.csv")

# data= mydata
data = data %>% select(-X)

# simulate some nco outcomes 
data$nco1 = rbinom(nrow(data),1,0.3)
data$nco2 = rbinom(nrow(data),1,0.3)
data$nco3 = rbinom(nrow(data), 1, 0.3)

# simulate the time to event for nco outcomes
data$nco1_time = sample(data$time, nrow(data), replace = TRUE)
data$nco2_time = sample(data$time, nrow(data), replace = TRUE)
data$nco3_time = sample(data$time, nrow(data), replace = TRUE)


LATTE_ADRD = data.frame(data)
# Separate data for LATTE analysis
site_data <- list(
    site1 = LATTE_ADRD[LATTE_ADRD$site == "site1", ],
    site2 = LATTE_ADRD[LATTE_ADRD$site == "site2", ],
    site3 = LATTE_ADRD[LATTE_ADRD$site == "site3", ],
    site4 = LATTE_ADRD[LATTE_ADRD$site == "site4", ],
    site5 = LATTE_ADRD[LATTE_ADRD$site == "site5", ],
    site6 = LATTE_ADRD[LATTE_ADRD$site == "site6", ],
    site7 = LATTE_ADRD[LATTE_ADRD$site == "site7", ],
    site8 = LATTE_ADRD[LATTE_ADRD$site == "site8", ],
    site9 = LATTE_ADRD[LATTE_ADRD$site == "site9", ],
    site10 = LATTE_ADRD[LATTE_ADRD$site == "site10", ]
)



# Non-Concurrent Outcomes (NCOs) setup for LATTE calibration
nco_outcomes = c(
    "nco1", "nco2", "nco3"
)

# --- Setup: Working directory ---
if (!dir.exists("pda_latte_results")) {
  dir.create("pda_latte_results")
}
original_wd <- getwd()
setwd("../pda/LATTE/")

# 1. Run Pooled Analysis (Traditional Benchmark)
cat("## Running Traditional Pooled Analysis...\n")
pooled_results <- run_pooled_analysis(LATTE_ADRD, outcome_id, outcome_time, sites, treatment_var, xvars)

# 2. Setup LATTE Control Parameters
### LATTE example, with NCO outcomes, logistic outcome model
control <- list(
    project_name = "LATTE Demo Study",
    step = "initialize", # Will be set by run_latte_analysis
    sites = sites,
    model = "LATTE",
    family = "binomial",
    outcome = outcome_id,
    nco_outcomes = nco_outcomes,
    variables = xvars,
    treatment = treatment_var,
    lead_site = "site1",
    balancing_method = "stratification", # stratification
    outcome_model = "logistic" # logistic
)

### LATTE example, without NCO outcomes 
# control <- list(
#     project_name = "LATTE Demo Study",
#     step = "initialize", # Will be set by run_latte_analysis
#     sites = sites,
#     model = "LATTE",
#     family = "binomial",
#     outcome = outcome_id,
#     variables = xvars,
#     treatment = treatment_var,
#     lead_site = "site1",
#     balancing_method = "stratification", # stratification
#     outcome_model = "logistic" # logistic
# )

### LATTE example, with NCO outcomes, poisson outcome model
# control <- list(
#     project_name = "LATTE Demo Study",
#     step = "initialize", # Will be set by run_latte_analysis
#     sites = sites,
#     model = "LATTE",
#     family = "binomial",
#     outcome = outcome_id,
#     nco_outcomes = nco_outcomes,
#     variables = xvars,
#     treatment = treatment_var,
#     lead_site = "site1",
#     balancing_method = "stratification", # stratification
#     outcome_model = "poisson", # logistic
#     outcome_times = c("time", "nco1_time", "nco2_time", "nco3_time")
# )


# 3. Run LATTE Analysis (Distributed Simulation)
cat("\n## Running Distributed LATTE Analysis (Initialize & Estimate)...\n")
dir.create("pda_latte_results")
setwd("pda_latte_results")
# Step 1: Initialize at each site (each site computes their 2x2 tables)
menu <- function(choices, title = NULL) 1

pda(site_id = 'site1', control = control, ipdata = site_data$site1, dir = getwd())
library(cobalt)
library(geex)
library(glmnet)
menu <- function(choices, title = NULL) 1

pda(site_id = 'site2', ipdata = site_data$site2, dir = getwd())
menu <- function(choices, title = NULL) 1

pda(site_id = 'site3', ipdata = site_data$site3, dir = getwd())
menu <- function(choices, title = NULL) 1

pda(site_id = 'site4', ipdata = site_data$site4, dir = getwd())
menu <- function(choices, title = NULL) 1
pda(site_id = "site5", ipdata = site_data$site5, dir = getwd())

menu <- function(choices, title = NULL) 1

pda(site_id = "site6", ipdata = site_data$site6, dir = getwd())
menu <- function(choices, title = NULL) 1
pda(site_id = "site7", ipdata = site_data$site7, dir = getwd())

menu <- function(choices, title = NULL) 1
pda(site_id = "site8", ipdata = site_data$site8, dir = getwd())
menu <- function(choices, title = NULL) 1
pda(site_id = "site9", ipdata = site_data$site9, dir = getwd())
menu <- function(choices, title = NULL) 1
pda(site_id = "site10", ipdata = site_data$site10, dir = getwd())


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


