require(survival)
require(data.table)
require(pda)
# data(lung)

## In the toy example below we aim to analyze the association of lung status with age and sex using Cox regression,
## data(lung) from 'survival', we randomly assign to 3 sites: 'site1', 'site2', 'site3'
## we demonstrate using PDA ODACH can obtain a surrogate estimator that is close to the pooled estimate.
## Different from ODAC, ODACH accounts for heterogeneity across sites by allowing site-specific baseline hazard functions and feature distributions.
## We run the example in local directory. In actual collaboration, account/password for pda server
## will be assigned to the sites at the server https://pda.one.
## Each site can access via web browser to check the communication of the summary stats.

data(lung2)
lung_split <- split(lung2, lung2$site)
## fit Cox PH reg using pooled data, stratified by site
fit.pool <- coxph(Surv(time, status) ~ age + sex + strata(site), data = lung2)

sites = c('site1', 'site2', 'site3')
S=readline(prompt="Type  <Return>   to continue : ")
# ############################  STEP 1: initialize  ###############################
control <- list(project_name = 'Lung cancer study',
                step = 'initialize',
                sites = sites,
                heterogeneity = TRUE,
                model = 'ODAC',
                family = 'cox',
                outcome = "Surv(time, status)",
                variables = c('age', 'sex'),
                # xlev = list(sex=c('F', 'M')),  #levels of all categorical X's, with the first being the reference
                optim_maxit = 100,
                lead_site = 'site1',
                upload_date = as.character(Sys.time()) )


## run the example in local directory:
## specify your working directory, default is the tempdir
mydir <- getwd()   # tempdir()
## assume lead site1: enter "1" to allow transferring the control file
pda(site_id = 'site1', control = control, dir = mydir)
## in actual collaboration, account/password for pda server will be assigned, thus:
# pda(site_id = 'site1', control = control, uri = 'https://pda.one', secret='abc123')
## you can also set your environment variables, and no need to specify them in pda:
# Sys.setenv(PDA_USER = 'site1', PDA_SECRET = 'abc123', PDA_URI = 'https://pda.one')
# pda(site_id = 'site1', control = control)

S=readline(prompt="Type  <Return>   to continue : ")
## assume remote site3: enter "1" to allow tranferring your local estimate
pda(site_id = 'site3', ipdata = lung_split[[3]], dir=mydir)

S=readline(prompt="Type  <Return>   to continue : ")
## assume remote site2: enter "1" to allow tranferring your local estimate
pda(site_id = 'site2', ipdata = lung_split[[2]], dir=mydir)

S=readline(prompt="Type  <Return>   to continue : ")
## assume lead site1: enter "1" to allow tranferring your local estimate
## control.json is also automatically updated
pda(site_id = 'site1', ipdata = lung_split[[1]], dir=mydir)

## if lead site1 initialized before other sites,
## lead site1: uncomment to sync the control before STEP 2
# pda(site_id = 'site1', control = control)
# config <- getCloudConfig(site_id = 'site1', dir=mydir)
# pdaSync(config)

S=readline(prompt="Type  <Return>   to continue : ")
# ###################   STEP 2: derivative  #################
## assume remote site3: enter "1" to allow tranferring your derivatives
pda(site_id = 'site3', ipdata = lung_split[[3]], dir=mydir)

S=readline(prompt="Type  <Return>   to continue : ")
## assume remote site2: enter "1" to allow tranferring your derivatives
pda(site_id = 'site2', ipdata = lung_split[[2]], dir=mydir)

S=readline(prompt="Type  <Return>   to continue : ")
## assume lead site1: enter "1" to allow tranferring your derivatives
pda(site_id = 'site1', ipdata = lung_split[[1]], dir=mydir)


S=readline(prompt="Type  <Return>   to continue : ")
# ############################  STEP 3: estimate  ###############################
## assume lead site1: enter "1" to allow tranferring the surrogate estimate
pda(site_id = 'site1', ipdata = lung_split[[1]], dir=mydir)

## the PDA ODACH is now completed!
## All the sites can still run their own surrogate estimates and broadcast them.

S=readline(prompt="Type  <Return>   to continue : ")

## compare the surrogate estimate with the pooled and meta estimates
config <- getCloudConfig(site_id = 'site1', dir=mydir)
fit.pda <- pdaGet(name = 'site1_estimate', config = config)
control <- pdaGet('control', config)
cbind(b.pool=fit.pool$coef,
      b.meta =control$beta_init,
      b.pda=fit.pda$btilde )

# S=readline(prompt="Type  <Return>   to continue : ")
# ## assume remote site2: (optional)
# pda(site_id = 'site2', ipdata = lung_split[[2]], dir=mydir)
# 
# S=readline(prompt="Type  <Return>   to continue : ")
# ## assume remote site3: (optional)
# pda(site_id = 'site3', ipdata = lung_split[[3]], dir=mydir)
# 
# 
# S=readline(prompt="Type  <Return>   to continue : ")
# ## If all the sites broadcast their surrogate estimates,
# ## a final synthesize step can further improve the estimate.
# ## assume lead site1: uncomment to synchoronize the control before STEP 4
# pda(site_id = 'site1', control = control, dir = mydir)
# config <- getCloudConfig(site_id = 'site1', dir = mydir)
# pdaSync(config)
# 
# S=readline(prompt="Type  <Return>   to continue : ")
# # ########################  STEP 4: synthesize (optional)  ########################
# ## assume lead site1:
# pda(site_id = 'site1', ipdata = lung_split[[1]], dir=mydir)
