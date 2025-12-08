require(survival)
require(data.table)
require(pda)

## In the toy example below we aim to analyze the association of lung cancer 
## with age and sex using Cox regression with time-varying effects,
## data(lung) from 'survival', we randomly assign to 3 sites: 'site1', 'site2', 'site3'
## we demonstrate using PDA ODACT can obtain a surrogate estimator that is close to the pooled estimate.
## We run the example in local directory. 
## In actual collaboration, the data communication can be done via the PDA_OTA platform https://pda-ota.pdamethods.org/
## Each site can access via web browser to transfer aggregate data and check the progress of the project.

# load('pda/data/lung2.rda') 
data(lung2)
lung_split <- split(lung2, lung2$site)
h = 200   # window bandwidth  
evalt = seq(0,by=200,len=5) # time points to estimate beta(t)
px = 2    # two coef: age and sex

## fit Cox with beta(t) using pooled data, stratified by site
fit.pool <- mycoxph_bt(lung2[,-1], site=lung2$site, fn=llpl_st, times=evalt, h=h, betabar=c(0,0), hessian=T)
 


sites = c('site1', 'site2', 'site3')
S=readline(prompt="Type  <Return>   to continue : ")
# ############################  STEP 1: initialize  ###############################
control <- list(project_name = 'Lung cancer study',
                step = 'initialize',
                sites = sites,
                heterogeneity = TRUE,
                model = 'ODACT',
                family = 'cox',
                outcome = "Surv(time, status)",
                variables = c('age', 'sex'),
                times = seq(0,by=200,len=5), # time points to estimate beta(t)  
                bandwidth = 200,
                optim_maxit = 100,
                lead_site = 'site1',
                init_method = "meta",
                upload_date = as.character(Sys.time()) )


## run the example in local directory:
## specify your working directory, default is the tempdir
mydir <- getwd()   # tempdir()
## assume lead site1: enter "1" to allow transferring the control file
pda(site_id = 'site1', control = control, dir = mydir)


S=readline(prompt="Type  <Return>   to continue : ")
## assume remote site3: enter "1" to allow transferring your local estimate
pda(site_id = 'site3', ipdata = lung_split[[3]], dir=mydir)

S=readline(prompt="Type  <Return>   to continue : ")
## assume remote site2: enter "1" to allow transferring your local estimate
pda(site_id = 'site2', ipdata = lung_split[[2]], dir=mydir)

S=readline(prompt="Type  <Return>   to continue : ")
## assume lead site1: enter "1" to allow transferring your local estimate
## control.json is also automatically updated
pda(site_id = 'site1', ipdata = lung_split[[1]], dir=mydir)


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

## the PDA ODACT is now completed! 

## compare the surrogate estimate with the pooled and meta estimates
config <- getCloudConfig(site_id = 'site1', dir=mydir)
fit.pda <- pdaGet(name = 'site1_estimate', config = config)
control <- pdaGet('control', config)
cbind(b.pool=fit.pool$b.pool[1,], # 1st coef: the age effect
      b.meta =control$beta_init[1,],
      b.pda=fit.pda$btilde[1,])
 
