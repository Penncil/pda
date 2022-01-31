require(data.table)
require(pda)
require(glmnet)
library(MASS) 
data(ODACAT_ordinal) # load simulated dataset
# load('data/ODACAT_ordinal.rda')

## In the toy example below we aim to analyze the association of a 3-category outcome 'y' 
## with 3 covariates 'X1', 'X2', 'X3'  using multinom-logistic (mlogit, if y is nominal) 
## or proportional odds logistic regression (polr, if y is ordinal).
## In this simulated data, we have 3 sites: 'site1', 'site2', 'site3'.
## We demonstrate using PDA ODACAT can obtain a surrogate estimator that is close to the pooled estimate. 
## We run the example in local directory. In actual collaboration, account/password for pda server 
## will be assigned to the sites at the server https://pda.one.
## Each site can access via web browser to check the communication of the summary stats.

sites <- unique(ODACAT_ordinal$id.site)
dd <- split(ODACAT_ordinal, ODACAT_ordinal$id.site)


S=readline(prompt="Type  <Return>   to continue : ")
# ############################  STEP 1: initialize  ###############################
control <- list(project_name = 'ODACAT Example',
                step = 'initialize',
                sites = sites,
                heterogeneity = FALSE,
                model = 'ODACAT',
                family = 'multicategory',
                number_outcome_categories=3,
                ordinal_categories=TRUE,
                outcome = "outcome",
                variables = c('X1','X2','X3'),
                # xlev = list(sex=c('F', 'M')),  #levels of all categorical X's, with the first being the reference
                optim_maxit = 100,
                lead_site = 'site1',
                upload_date = as.character(Sys.time()) )

## fit polr reg using pooled data
ODACAT_ordinal$outcome=factor(ODACAT_ordinal$outcome,levels=c(1,2,3))
fit.pool.polr=polr(outcome~X1+X2+X3,data=ODACAT_ordinal)

## run the example in local directory:
## specify your working directory, default is the tempdir
mydir <- getwd()  # tempdir()
## assume lead site1: enter "1" to allow transferring the control file 
pda(site_id = 'site1', control = control, dir = mydir)


## in actual collaboration, account/password for pda server will be assigned, thus:
# pda(site_id = 'site1', control = control, uri = 'https://pda.one', secret='abc123')
## you can also set your environment variables, and no need to specify them in pda:
# Sys.setenv(PDA_USER = 'site1', PDA_SECRET = 'abc123', PDA_URI = 'https://pda.one')
# pda(site_id = 'site1', control = control)

S=readline(prompt="Type  <Return>   to continue : ")
## assume remote site3: enter "1" to allow tranferring your local estimate 
pda(site_id = 'site3', ipdata = dd[[3]][,-1], dir=mydir)

S=readline(prompt="Type  <Return>   to continue : ")
## assume remote site2: enter "1" to allow tranferring your local estimate  
pda(site_id = 'site2', ipdata = dd[[2]][,-1], dir=mydir)

S=readline(prompt="Type  <Return>   to continue : ")
## assume lead site1: enter "1" to allow tranferring your local estimate  
## control.json is also automatically updated
pda(site_id = 'site1', ipdata = dd[[2]][,-1], dir=mydir)

## if lead site1 initialized before other sites,
## lead site1: uncomment to sync the control before STEP 2
# pda(site_id = 'site1', control = control)
# config <- getCloudConfig(site_id = 'site1')
# pdaSync(config)

S=readline(prompt="Type  <Return>   to continue : ")
# ############################  STEP 2: derivative  ###############################
## assume remote site3: enter "1" to allow tranferring your derivatives  
pda(site_id = 'site3', ipdata = dd[[3]][,-1], dir=mydir)

S=readline(prompt="Type  <Return>   to continue : ")
## assume remote site2: enter "1" to allow tranferring your derivatives  
pda(site_id = 'site2', ipdata = dd[[2]][,-1], dir=mydir)

S=readline(prompt="Type  <Return>   to continue : ")
## assume lead site1: enter "1" to allow tranferring your derivatives  
pda(site_id = 'site1', ipdata = dd[[1]][,-1], dir=mydir)

S=readline(prompt="Type  <Return>   to continue : ")
# ############################  STEP 3: estimate  ###############################
## assume lead site1: enter "1" to allow tranferring the surrogate estimate  
pda(site_id = 'site1', ipdata = dd[[1]][,-1], dir=mydir)

## the PDA ODACAT is now completed!

S=readline(prompt="Type  <Return>   to continue : ")


## compare the surrogate estimate with the pooled estimate
config <- getCloudConfig(site_id = 'site1', dir=mydir)
fit.odact <- pdaGet(name = 'site1_estimate', config = config)
cbind(b.pool=c(fit.pool.polr$coef, fit.pool.polr$zeta), b.odact=fit.odact$btilde)

# # res.meta=pdaGet(name = 'control', config = config) 
# res.meta <- meta.fit(p=3,control,config=config)
# tt = repar(res.meta$bbar, diag(res.meta$vmeta), 3)
# cbind(b.pool=c(fit.pool.polr$coef, fit.pool.polr$zeta),
#       b.meta=c(tt$beta, tt$zeta), 
#       b.odact=fit.odact$btilde)

