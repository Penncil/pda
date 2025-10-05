require(lme4)
require(minqa)
require(data.table)
require(geex)
## In the toy example below we aim to estimate the average treatment effect (ATE) of covid vaccination 
# on the outcome (number of Post Acute Sequelae of COVID features), with 50 covariates X1-X50
# The data are simulated as from 3 sites (site1, site2, site3), each has 300 participants.

# library(devtools)
# devtools::install_github("penncil/pda")
# library(pda)
library(dplyr)
library(numDeriv)
library(glmnet)

 
Rcpp::sourceCpp("pda/src/DisC2o.cpp")
# load("~/Dropbox/PDA-git/pda/data/long_covid.rda")
data(long_covid)   
long_covid_split <- split(long_covid, long_covid$site)
 

#############################################
############ PS model fitting ###############
#############################################

# ############################  STEP 1: initialize  ###############################
variables <- names(DisC2o_data)[4:53]
sites <- unique(DisC2o_data$site)
control <- list(project_name = 'PASC_vaccination',
                sites = sites,
                lead_site = 'site1', 
                step = 'PSinitialize',
                init_method = "lead",
                heterogeneity = FALSE,
                model = 'DisC2o',
                family = 'guassian',
                treatment='covid_vaccination', 
                outcome = 'PASC_features', 
                variables = variables,
                optim_maxit = 100,
                upload_date = as.character(Sys.time()) )


## run the example in local directory:
## specify your working directory, default is the tempdir
mydir <- getwd()  # tempdir()
## assume lead site1: enter "1" to allow transferring the control file 
## PDA-OTA interactive section: use pda-ota for data communication
# STEP 0 [ALL sites]: remove any json files if exist
file.remove(list.files(mydir,full.names = T)[grepl('.json', list.files(mydir))])

# STEP 0 [LEAD site only]: create control.json, read the R output
pda(site_id = control$lead_site, control=control, dir=mydir)
# and upload control.json to pda-ota 


S=readline(prompt="Type  <Return>   to continue : ")
# ############################ PS model estimate  ###############################
# Step 1 [LEAD site only]: 
# calculate (individual PS model estimates), read the R output   
pda(site_id = control$lead_site, ipdata = long_covid_split[[1]], dir=mydir)
# and upload *_initialize.json to pda-ota


S=readline(prompt="Type  <Return>   to continue : ")
# STEP 2 [ALL sites]: download control.json to your working dir,
# calculate PS model first and second order derivatives of other sites, 
pda(site_id = 'site3', ipdata = long_covid_split[[3]], dir=mydir)


S=readline(prompt="Type  <Return>   to continue : ")
# STEP 2 [ALL sites]: download control.json to your working dir,
# calculate PS model first and second order derivatives of other sites,
pda(site_id = 'site2', ipdata = long_covid_split[[2]], dir=mydir)


S=readline(prompt="Type  <Return>   to continue : ")
# STEP 3 [LEAD site only]: 
# 1. download all (*_derive.json) to your working dir,
# 2. calculate PS model first and second order derivatives of lead site
# 3. calculate PS model estimate
# 4. calculate initial estimate for OM model
pda(site_id = control$lead_site, ipdata = long_covid_split[[1]], dir=mydir)
# close the project and upload the final result *_estimate.json 



# ############################   OM model estimate  ###############################
S=readline(prompt="Type  <Return>   to continue : ")
# STEP 4 [ALL sites]: download control.json to your working dir,
# calculate OM model first and second order derivatives of other sites,   
pda(site_id = 'site3', ipdata = long_covid_split[[3]], dir=mydir)

S=readline(prompt="Type  <Return>   to continue : ")
# STEP 4 [ALL sites]: download control.json to your working dir,
# calculate OM model first and second order derivatives of other sites, 
pda(site_id = 'site2', ipdata = long_covid_split[[2]], dir=mydir)


S=readline(prompt="Type  <Return>   to continue : ")
# STEP 5 [LEAD site only]: 
# 1. download all (*_derive.json) to your working dir,
# 2. calculate PS model first and second order derivatives of lead site
# 3. calculate PS model estimate
pda(site_id = control$lead_site, ipdata = long_covid_split[[1]], dir=mydir)
# close the project and upload the final result *_estimate.json 

# ############################   AIPW  estimate  ###############################
S=readline(prompt="Type  <Return>   to continue : ")
# STEP 6 [all site]: download all files (*.json) to your working dir,
# calculate AIPW estimate for each site  
pda(site_id = 'site3', ipdata = long_covid_split[[3]], dir=mydir)

S=readline(prompt="Type  <Return>   to continue : ")
# STEP 6 [all site]: download all files (*.json) to your working dir,
# calculate AIPW estimate for each site 
pda(site_id = 'site2', ipdata = long_covid_split[[2]], dir=mydir)


S=readline(prompt="Type  <Return>   to continue : ")
# STEP 10 [LEAD site only]: 
# 1. calculate AIPW estimate for lead site 
# 2. calculate aggregated AIPW estimate
pda(site_id = control$lead_site, ipdata = long_covid_split[[1]], dir=mydir)
# close the project and upload the final result *_estimate.json 


