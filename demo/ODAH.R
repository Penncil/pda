# require(countreg)
require(data.table)
require(pda)
 
## In the toy example below we aim to analyze the association of CrabSatellites with width and weight using hurdle regression,
## data(CrabSatellites) from 'countreg', we randomly assign to 2 sites: 'site1', 'site2' 
## we demonstrate using PDA ODAH, can obtain a surrogate estimator that is close to the pooled estimate. 
## We run the example in local directory. In actual collaboration, account/password for pda server 
## will be assigned to the sites at the server https://pda.one.
## Each site can access via web browser to check the communication of the summary stats.


# install.packages("countreg", repos="http://R-Forge.R-project.org")
require(countreg)
# Create 2 sites, split the  data amongst them
sites = c('site1', 'site2')
set.seed(42)
data("CrabSatellites", package = "countreg")
cs <- CrabSatellites[, c("satellites", "width", "weight")]
dd_split <- split(cs, sample(1:length(sites), nrow(cs), replace=TRUE))

## fit logistic reg using pooled data
fit.pool <- hurdle(satellites ~ width+weight, data = cs)

 
S=readline(prompt="Type  <Return>   to continue : ")
# ############################  STEP 1: initialize  ###############################
control <- list(project_name = 'CrabSatellites study',
                step = 'initialize',
                sites = sites,
                heterogeneity = FALSE,
                model = 'ODAH',
                family = 'hurdle',
                outcome = "satellites",
                variables_hurdle_count = c("width", "weight"),
                variables_hurdle_zero = c("width", "weight"),
                offset = NULL,
                optim_maxit = 100,
                lead_site = 'site1',
                upload_date = as.character(Sys.time()) )


## run the example in local directory:
## specify your working directory, default is the tempdir
mydir <- getwd()   #  tempdir()
## assume lead site1: enter "1" to allow transferring the control file  
pda(site_id = 'site1', control = control, dir = mydir)
## in actual collaboration, account/password for pda server will be assigned, thus:
# pda(site_id = 'site1', control = control, uri = 'https://pda.one', secret='abc123')
## you can also set your environment variables, and no need to specify them in pda:
# Sys.setenv(PDA_USER = 'site1', PDA_SECRET = 'abc123', PDA_URI = 'https://pda.one')
# pda(site_id = 'site1', control = control)

S=readline(prompt="Type  <Return>   to continue : ")
## assume remote site2: enter "1" to allow tranferring your local estimate  
pda(site_id = 'site2', ipdata = dd_split[[2]], dir=mydir)


S=readline(prompt="Type  <Return>   to continue : ")
## assume lead site1: enter "1" to allow tranferring your local estimate  
## control.json is also automatically updated
pda(site_id = 'site1', ipdata = dd_split[[1]], dir=mydir)

## if lead site1 initialized before other sites,
## lead site1: uncomment to sync the control before STEP 2
# pda(site_id = 'site1', control = control)
# config <- getCloudConfig(site_id = 'site1')
# pdaSync(config)

S=readline(prompt="Type  <Return>   to continue : ")
# ############################  STEP 2: derivative  ###############################
## assume remote site2: enter "1" to allow tranferring your derivatives  
pda(site_id = 'site2', ipdata = dd_split[[2]], dir=mydir)

S=readline(prompt="Type  <Return>   to continue : ")
## assume lead site1: enter "1" to allow tranferring your derivatives  
pda(site_id = 'site1', ipdata = dd_split[[1]], dir=mydir)

S=readline(prompt="Type  <Return>   to continue : ")
# ############################  STEP 3: estimate  ###############################
## assume lead site1: enter "1" to allow tranferring the surrogate estimate  
pda(site_id = 'site1', ipdata = dd_split[[1]], dir=mydir)

## the PDA ODAH is now completed!
## All the sites can still run their own surrogate estimates and broadcast them.

S=readline(prompt="Type  <Return>   to continue : ")
## compare the surrogate estimate with the pooled estimate
config <- getCloudConfig(site_id = 'site1', dir=mydir)
fit.odah <- pdaGet(name = 'site1_estimate', config = config)
control <- pdaGet(name = 'control', config = config)
cbind(b.count.pool=fit.pool$coef$count,
      b.count.meta=control$beta_count_init,
      b.count.odah=fit.odah$btilde_count,
      b.zero.pool=fit.pool$coef$zero,
      b.zero.meta=control$beta_zero_init,
      b.zero.odah=fit.odah$btilde_zero )

# S=readline(prompt="Type  <Return>   to continue : ")
# ## assume remote site2: (optional)
# pda(site_id = 'site2', ipdata = dd_split[[2]], dir=mydir)
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
# pda(site_id = 'site1', ipdata = dd_split[[1]], dir=mydir)
