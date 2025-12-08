require(lme4)
require(minqa)
require(data.table)
require(geex)
require(pda) 
library(dplyr)
library(numDeriv)

## In the toy example below we aim to analyze the association of COVID-19 status
## with age, sex and lab test using linear model, 
## data(sample_data_covid) is simulated and assumed to come from 3 sites: 'site1', 'site2', 'site3'
## we demonstrate using PDA COLA can obtain identical estimation as the pooled analysis.
## 3 models are tested:
# - COLA-GLM
# - COLA-GLM-H
# - COLA-GLMM
## We run the example in local directory. 
## In actual collaboration, the data communication can be done via the PDA_OTA platform https://pda-ota.pdamethods.org/
## Each site can access via web browser to transfer aggregate data and check the progress of the project.
# library(devtools)
# devtools::install_github("penncil/pda")
 
 
  
data("COLA_covid", package = "pda")

sites = unique(COLA_covid$site) # c('site1', 'site2', 'site3')
data_split <- split(COLA_covid, COLA_covid$site)

S=readline(prompt="Type  <Return>   to continue : ")
# # ########################  COLA-GLM ############################################
## binary outcome
# ############################  STEP 1: initialize  ###############################
control <- list(project_name = 'COVID-19 study',
                step = 'initialize',
                sites = sites,
                heterogeneity = FALSE,
                mixed_effects = FALSE,
                model = 'COLA',
                family = 'binomial',
                outcome = "status",
                variables = c('age', 'sex', 'medical_condition'),
                cutoff = NULL, # set this cutoff as NULL or set a number (e.g., 5 or 11)
                link = "canonical", 
                optim_maxit = 100,
                lead_site = 'site1',
                upload_date = as.character(Sys.time()))


## run the example in local directory:
## specify your working directory, default is the tempdir
mydir <- getwd()  # tempdir()
## assume lead site1: enter "1" to allow transferring the control file 
pda(site_id = 'site1', control = control, dir = mydir)
  

S=readline(prompt="Type  <Return>   to continue : ")
## assume remote site3: enter "1" to allow tranferring your local estimate 
pda(site_id = 'site3', ipdata = data_split[[3]], dir=mydir)
 
S=readline(prompt="Type  <Return>   to continue : ")
## assume remote site2: enter "1" to allow tranferring your local estimate  
pda(site_id = 'site2', ipdata = data_split[[2]], dir=mydir)
 
S=readline(prompt="Type  <Return>   to continue : ")
## assume lead site1: enter "1" to allow tranferring your local estimate  
## control.json is also automatically updated
pda(site_id = 'site1', ipdata = data_split[[1]], dir=mydir)


S=readline(prompt="Type  <Return>   to continue : ")
# ############################  STEP 2: estimate  ###############################
## assume lead site1: enter "1" to allow tranferring the surrogate estimate  
pda(site_id = 'site1', ipdata = data_split[[1]], dir=mydir)
 
## the PDA COLA is now completed!
