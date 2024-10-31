require(survival)
require(data.table)
# require(pda) 

# rho = 0.2 
# odach_cc$sampling_weight[synthetic_data$subcohort==1] = 1/rho  
  
# load('/Users/chongliangluo/Library/CloudStorage/Dropbox/PDA-git/pda/data/odach_cc.rda')
# setwd('/Users/chongliangluo/Library/CloudStorage/Dropbox/PDA_test/CL/')

## In the toy example below we aim to analyze the association of survival {time} with X1 and X2 using Cox reg with case-cohort design,
## data(odach_cc) simulated, subsampled case+cohort data, 3 sites: 'site1', 'site2', 'site3'
## the full_cohort_size are 800 600 400 for the 3 sites respectively 
## we demonstrate using PDA ODACH_CC can obtain a surrogate estimator that is close to the pooled estimate.
## ODACH_CC relies on the surrogate of the pooled stratified case-cohort pseudo likelihood function
## We run the example in local directory. In actual collaboration, account/password for pda server
## will be assigned to the sites at the server https://pda.one.
## Each site can access via web browser to check the communication of the summary stats.

data(odach_cc)
data_split <- split(odach_cc, odach_cc$site)

# ## cch stratified by site...
# fit.strat <- optim(par = initial_beta, fn = pool_fun, control = list(fnscale = -1), method = "BFGS",
#                 covariate_list = pre_processing$covariate_list,
#                 failure_position = pre_processing$failure_position,
#                 failure_num = pre_processing$failure_num,
#                 risk_sets = pre_processing$risk_sets,
#                 risk_set_weights = pre_processing$risk_set_weights,
#                 K = length(data_list))


## fit Cox case-cohort using survival::cch with pooled data, NO stratified by site
fit.pool <- cch(Surv(time, status) ~ X1 + X2, data = cbind(ID=1:nrow(odach_cc), odach_cc),
                subcoh = ~subcohort, id = ~ID, cohort.size = 800+600+400, method = 'Prentice')
# coxph(Surv(time, status) ~ X1 + X2+ strata(site), data = odach_cc )$coef
fit.pool$coef 
fit.pool$var     
    

# cch each site
cch(Surv(time, status) ~ X1 + X2, data = cbind(ID=1:nrow(data_split[[3]]), data_split[[3]]),
    subcoh = ~subcohort, id = ~ID, cohort.size = 400, method = 'Prentice')$coef
# -0.6578062  0.6004303 
# -0.5175115  0.3163663 
# -0.5471417  0.4696186 

sites = c('site1', 'site2', 'site3')
S=readline(prompt="Type  <Return>   to continue : ")
# ############################  STEP 1: initialize  ###############################
control <- list(project_name = 'ODACH case-cohort toy example',
                step = 'initialize',
                sites = sites,
                heterogeneity = TRUE,
                model = 'ODACH_CC',
                family = 'cox',
                full_cohort_size = c(site1=800,site2=600,site3=400),
                method = 'Prentice',
                outcome = "Surv(time, status)",
                variables = c('X1', 'X2'),
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
pda(site_id = 'site3', ipdata = data_split[[3]], dir=mydir) 
  
S=readline(prompt="Type  <Return>   to continue : ")
## assume remote site2: enter "1" to allow tranferring your local estimate
pda(site_id = 'site2', ipdata = data_split[[2]], dir=mydir)

S=readline(prompt="Type  <Return>   to continue : ")
## assume lead site1: enter "1" to allow tranferring your local estimate
## control.json is also automatically updated
pda(site_id = 'site1', ipdata = data_split[[1]], dir=mydir)
# "bbar":[-0.55271374157,0.41587314193]

## if lead site1 initialized before other sites,
## lead site1: uncomment to sync the control before STEP 2
# pda(site_id = 'site1', control = control)
# config <- getCloudConfig(site_id = 'site1', dir=mydir)
# pdaSync(config)

S=readline(prompt="Type  <Return>   to continue : ")
# ###################   STEP 2: derivative  #################
## assume remote site3: enter "1" to allow tranferring your derivatives
pda(site_id = 'site3', ipdata = data_split[[3]], dir=mydir)
# "logL_D1":[0.23547460153,0.92154395093],"logL_D2":[[-23.4031336829,-2.0585355377],[-2.0585355377,-17.1356300437]]

S=readline(prompt="Type  <Return>   to continue : ")
## assume remote site2: enter "1" to allow tranferring your derivatives
pda(site_id = 'site2', ipdata = data_split[[2]], dir=mydir)
# "logL_D1":[-0.68819757269,1.9822067467],"logL_D2":[[-16.342905835,-5.4532522136],[-5.4532522136,-14.2454398028]]

S=readline(prompt="Type  <Return>   to continue : ")
## assume lead site1: enter "1" to allow tranferring your derivatives
pda(site_id = 'site1', ipdata = data_split[[1]], dir=mydir)
# "logL_D1":[1.7694795086,-3.3851725372],"logL_D2":[[-41.3249413889,3.005647082],[3.005647082,-32.9944494016]]

S=readline(prompt="Type  <Return>   to continue : ")
# ############################  STEP 3: estimate  ###############################
## assume lead site1: enter "1" to allow tranferring the surrogate estimate
pda(site_id = 'site1', ipdata = data_split[[1]], dir=mydir)
 

## the PDA ODACH_CC is now completed!
## All the sites can still run their own surrogate estimates and broadcast them.
 
