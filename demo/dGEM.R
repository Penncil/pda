require(survival)
require(data.table)
require(pda)


## In the toy example below we aim to conduct hospital profiling using dGEM. We randomly assign to 3 sites: 'site1', 'site2', 'site3'
## we demonstrate using PDA ODACH can obtain a surrogate estimator that is close to the pooled estimate.
## Different from ODAC, ODACH accounts for heterogeneity across sites by allowing site-specific baseline hazard functions and feature distributions.
## We run the example in local directory. In actual collaboration, account/password for pda server
## will be assigned to the sites at the server https://pda.one.
## Each site can access via web browser to check the communication of the summary stats.

data(lung2)
lung_split <- split(lung2, lung2$site)

sites = c('site1', 'site2', 'site3')
S=readline(prompt="Type  <Return>   to continue : ")
# ############################  STEP 1: initialize  ###############################
control <- list(project_name = 'Lung cancer study',
                step = 'initialize',
                sites = c('site1','site2','site3'),
                heterogeneity = TRUE,
                model = 'dGEM',
                family = 'binomial',
                outcome = "status",
                variables = c('age', 'sex'),
                variables_site_level = c('volume'),
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
# ###################   STEP 2: derive  #################
## assume remote site3: enter "1" to allow tranferring your derivatives
##' run dGEM step 2 under site 3
pda(site_id = 'site3', ipdata = lung_split[[3]], hosdata = c('volume' = 300), dir=getwd())

S=readline(prompt="Type  <Return>   to continue : ")
## assume remote site2: enter "1" to allow tranferring your derivatives
##' run dGEM step 2 under site 2
pda(site_id = 'site2', ipdata = lung_split[[2]], hosdata = c('volume' = 270), dir=getwd())

S=readline(prompt="Type  <Return>   to continue : ")
## assume lead site1: enter "1" to allow tranferring your derivatives
##' run dGEM step 2 under site 1
##' control.json is also automatically updated
pda(site_id = 'site1', ipdata = lung_split[[1]], hosdata = c('volume' = 150), dir=getwd())


S=readline(prompt="Type  <Return>   to continue : ")
# ############################  STEP 3: estimate  ###############################
## assume lead site1: enter "1" to allow tranferring the surrogate estimate
##' run dGEM step 3 under site 3
pda(site_id = 'site3', ipdata = lung_split[[3]], dir=getwd())

S=readline(prompt="Type  <Return>   to continue : ")
## assume lead site1: enter "1" to allow tranferring the surrogate estimate
##' run dGEM step 3 under site 2
pda(site_id = 'site2', ipdata = lung_split[[2]], dir=getwd())

S=readline(prompt="Type  <Return>   to continue : ")
## assume lead site1: enter "1" to allow tranferring the surrogate estimate
##' run dGEM step 3 under site 1
##' control.json is also automatically updated
pda(site_id = 'site1', ipdata = lung_split[[1]], dir=getwd())



S=readline(prompt="Type  <Return>   to continue : ")
#' ############################'  STEP 4: synthesize  ###############################
##' assume lead site1: enter "1" to allow tranferring the surrogate estimate
##' # step for lead site only. Other sites don't need to do this
pda(site_id = 'site1', ipdata = lung_split[[1]], dir=getwd())


S=readline(prompt="Type  <Return>   to continue : ")
# FINAL results
config <- getCloudConfig(site_id = 'site1', dir=getwd())
dGEM_event_rate <- pdaGet(name = 'site1_synthesize', config = config)
dGEM_event_rate
##' the PDA dGEM is now completed!
