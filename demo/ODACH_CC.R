require(survival)
require(data.table)
# require(pda) 
require(spatstat) # for using weighted.median() as init

# load('/Users/chongliang/Dropbox/PDA-git/pda/data/odach_cc.rda')
# setwd('/Users/chongliang/Dropbox/PDA_test/CL/')

## In the toy example below we aim to analyze the association of survival {time} 
##  with X1, Category, and Group using Cox reg with case-cohort design,
## data(odach_cc) simulated, subsampled case+cohort data, 3 sites: 'site1', 'site2', 'site3'
## the full_cohort_size are 800 600 400 for the 3 sites respectively 
## we demonstrate using PDA ODACH_CC can obtain a surrogate estimator that is close to the pooled estimate.
## ODACH_CC relies on the surrogate of the pooled stratified case-cohort pseudo likelihood function
## We run the example in local directory. In actual collaboration, account/password for pda server
## will be assigned to the sites at the server https://pda.one.
## Each site can access via web browser to check the communication of the summary stats.

# data(odach_cc)
odach_cc$Category = as.character(odach_cc$Category)
# to test site3 has missing category Z for "Category" variable
odach_cc[odach_cc$site=='site3'&odach_cc$Category=='Z (X,Y,Z)', ]$Category='Y (X,Y,Z)'
table(odach_cc$site, odach_cc$Category)
data_split <- split(odach_cc, odach_cc$site)

sites = c('site1', 'site2', 'site3')
control <- list(project_name = 'ODACH case-cohort toy example',
                step = 'initialize',
                sites = sites,
                heterogeneity = TRUE,
                model = 'ODACH_CC',
                family = 'cox',
                outcome = "Surv(time, status)",
                variables = c('X1','`Group (A,B,C)`', 'Category'),
                #levels of all categorical X's, with the first being the reference
                variables_lev = list(`Group (A,B,C)`=c('A','B','C'), 
                                     Category=c('X (X,Y,Z)','Y (X,Y,Z)','Z (X,Y,Z)')),  
                full_cohort_size = c(site1=800,site2=600,site3=400), # for ODACH_CC
                method = 'Prentice',                          # for ODACH_CC weights
                init_method = 'meta', # 'meta','median','weighted.median','lead'
                optim_method = 'BFGS',
                optim_maxit = 100,
                lead_site = 'site1',
                upload_date = as.character(Sys.time()) )

## cch stratified by site, using self-written cch_pooled()
fit.pool <- cch_pooled(Surv(time, status) ~ X1+`Group (A,B,C)`+Category, data=odach_cc, variables_lev = control$variables_lev,
                       full_cohort_size=c(site1=800,site2=600,site3=400), method='Prentice',var_sandwich=T)
fit.pool$par
#          X1   `Group (A,B,C)`B  `Group (A,B,C)`C CategoryY (X,Y,Z) CategoryZ (X,Y,Z) 
# -0.5129684         0.1277404         0.2421812        -0.2374705         0.1359541 
sqrt(diag(fit.pool$var)) # robust s.e. from sandwich 
#  0.1153119         0.3363288         0.3260169         0.3059613         0.3590618  

## another approach to do pooled: Vivian's coxph trick
precision <- min(diff(sort(odach_cc$time)))/2 # 1e-4
odach_cc$time_in = 0
odach_cc[odach_cc$subcohort == 0, "time_in"] <- odach_cc[odach_cc$subcohort == 0, "time"] - precision
odach_cc$ID = 1:nrow(odach_cc)
fit.pool1 <- coxph(Surv(time_in, time, status) ~ X1+`Group (A,B,C)`+Category + strata(site)+cluster(ID), data=odach_cc, robust=T)
fit.pool1$coef
# -0.5129683         0.1277399         0.2421812        -0.2374701         0.1359542
sqrt(diag(fit.pool1$var)) # identical to the above cch_pooled
# 0.1153119 0.3363288 0.3260168 0.3059612 0.3590618
# summary(fit.pool1)$coef[,'se(coef)']
# 0.1105210         0.3184666         0.2969596         0.2831813         0.3211316


# ############################  STEP 1: initialize  ############################### 
control$init_method = 'meta' # 'lead' # 'weighted.median' #  'median'

## run the example in local directory:
## specify your working directory, default is the tempdir
mydir <- getwd()   # tempdir()
file.remove(list.files(mydir,full.names = T)[grepl('.json', list.files(mydir))])

pda(site_id = 'site1', control = control, dir = mydir, upload_without_confirm = T,silent_message =F)

pda(site_id = 'site3', ipdata = data_split[[3]], dir=mydir, upload_without_confirm = T)
  
pda(site_id = 'site2', ipdata = data_split[[2]], dir=mydir, upload_without_confirm = T)

## control.json is also automatically updated
pda(site_id = 'site1', ipdata = data_split[[1]], dir=mydir, upload_without_confirm = T)

# ###################   STEP 2: derivative  #################
pda(site_id = 'site3', ipdata = data_split[[3]], dir=mydir,upload_without_confirm=T)

pda(site_id = 'site2', ipdata = data_split[[2]], dir=mydir,upload_without_confirm=T)

pda(site_id = 'site1', ipdata = data_split[[1]], dir=mydir,upload_without_confirm=T)

# ############################  STEP 3: estimate  ###############################
pda(site_id = 'site1', ipdata = data_split[[1]], dir=mydir,upload_without_confirm=T)
# "beta_init":[-0.5169,0.1755,0.2722,-0.2516,-0.0453] # (meta)
# "btilde":   [-0.5131,0.1275,0.2419,-0.2377,0.1373],
# "setilde":  [ 0.1156,0.3363,0.3312, 0.2970,0.3688]
              

############ 20250311:  to get better robust s.e. estimate ###################### 
# run one more round of ODACH_CC use the obtained ODACH_CC point est as initial
# reset control step and beta_init
config <- getCloudConfig(site_id =control$lead_site, dir=mydir)
fit.pda = pdaGet(name = 'site1_estimate', config = config)
control = pdaGet(name = 'control', config = config)
control$step = 'derive'
control$beta_init = fit.pda$btilde
pdaPut(control,'control', config = config, upload_without_confirm = T, silent_message = F)

# ###################   STEP 4: derivative again #################
pda(site_id = 'site3', ipdata = data_split[[3]], dir=mydir,upload_without_confirm=T)

pda(site_id = 'site2', ipdata = data_split[[2]], dir=mydir,upload_without_confirm=T)

pda(site_id = 'site1', ipdata = data_split[[1]], dir=mydir,upload_without_confirm=T)

# ############################  STEP 5: estimate again  ###############################
pda(site_id = 'site1', ipdata = data_split[[1]], dir=mydir,upload_without_confirm=T)
# "btilde": [-0.5086,0.1178,0.2315,-0.2312,0.2533]
# "setilde":[ 0.1155,0.3364,0.3278, 0.3004,0.3645]
             
fit.pda = pdaGet(name = 'site1_estimate', config = config)

## all coef est:
# pooled:                     [-0.5130,0.1277,0.2422,-0.2375, 0.1360] 
# meta:                       [-0.5169,0.1755,0.2722,-0.2516,-0.0453]
# odach_cc(meta):             [-0.5131,0.1275,0.2419,-0.2377, 0.1373] # first round 
# odach_cc(again):            [-0.5086,0.1178,0.2315,-0.2312, 0.2533] # second round 

## all s.e. est: 
# cch_pooled  0.1153 0.3363 0.3260 0.3060 0.3591      
# "setilde": [0.1156,0.3363,0.3312, 0.2970,0.3688]  # first round
# "setilde": [0.1155,0.3364,0.3278, 0.3004,0.3645]    # second round
# btilde and setilde (robust) are almost identical to pooled cch! 


# if you want to check ODACH_CC s.e. from inv-Hessian...
# config <- getCloudConfig(site_id =control$lead_site, dir=mydir)
# fit.pda = pdaGet(name = 'site1_estimate', config = config)
# sqrt(diag(solve(fit.pda$Htilde)))  


## the PDA ODACH_CC is now completed!
#NOT TESTED!# All the sites can still run their own surrogate estimates and broadcast them.
 













# cch each site
sapply(1:3, function(i) cch(Surv(time, status) ~ X1+`Group (A,B,C)`+Category, 
                            data = cbind(ID=1:nrow(data_split[[i]]), data_split[[i]]),
                            subcoh = ~subcohort, id = ~ID, 
                            cohort.size =c(site1=800,site2=600,site3=400)[i], 
                            method = 'Prentice')$coef) 
#        XX1  X`Group (A,B,C)`B  X`Group (A,B,C)`C XCategoryY (X,Y,Z) XCategoryZ (X,Y,Z) 
# -0.5184543          0.4805595          0.8929081         -0.5006887         -0.1715390 
# -0.4707972          0.2383644          0.6133549         -0.7252045          0.1941511 
# -0.5406719         -0.1171266         -1.2518623          0.3119451 
