require(survival)
require(data.table)
# require(pda) 
require(spatstat) # for using weighted.median() as init

# load('/Users/chongliangluo/Dropbox/PDA-git/pda/data/odach_cc.rda')
# setwd('/Users/chongliangluo/Dropbox/PDA_test/CL/')

# n_rows <- nrow(odach_cc)
# odach_cc$Category <- factor(rep(c("X (X,Y,Z)", "Y (X,Y,Z)", "Z (X,Y,Z)"), length.out = n_rows))
# names(odach_cc)[8] = 'Group (A,B,C)'
# save(odach_cc, file='/Users/chongliangluo/Library/CloudStorage/Dropbox/PDA-git/pda/data/odach_cc.rda')


## In the toy example below we aim to analyze the association of survival {time} 
##  with X1, Category, and Group using Cox reg with case-cohort design,
## data(odach_cc) simulated, subsampled case+cohort data, 3 sites: 'site1', 'site2', 'site3'
## the full_cohort_size are 800 600 400 for the 3 sites respectively 
## we demonstrate using PDA ODACH_CC can obtain a surrogate estimator that is close to the pooled estimate.
## ODACH_CC relies on the surrogate of the pooled stratified case-cohort pseudo likelihood function
## We run the example in local directory. In actual collaboration, account/password for pda server
## will be assigned to the sites at the server https://pda.one.
## Each site can access via web browser to check the communication of the summary stats.

data(odach_cc)
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
sqrt(diag(-solve(fit.pool$hessian)))  # naive s.e. from inv hessian 
# 0.1105210         0.3184666         0.2969597         0.2831813         0.3211316 
sqrt(diag(fit.pool$var)) # robust s.e. from sandwich 
# 0.1358760         0.3567799         0.3109574         0.3122135         0.3692968 
# 0.1150233         0.3576904         0.3655076         0.3193014         0.3499519  # updated...

# precision <- min(diff(sort(odach_cc$time)))/2 # 1e-4
# odach_cc$time_in = 0
# odach_cc[odach_cc$subcohort == 0, "time_in"] <- odach_cc[odach_cc$subcohort == 0, "time"] - precision
# odach_cc$ID = 1:nrow(odach_cc)
# fit.pool <- coxph(Surv(time_in, time, status) ~ X1+`Group (A,B,C)`+Category + strata(site)+cluster(ID), data=odach_cc, robust=T)
# sqrt(diag(fit.pool$var))
#  0.1153119 0.3363288 0.3260168 0.3059612 0.3590618


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


# S=readline(prompt="Type  <Return>   to continue : ")
# ############################  STEP 1: initialize  ############################### 
control$init_method = 'meta' # 'lead' # 'weighted.median' #  'median'

## run the example in local directory:
## specify your working directory, default is the tempdir
mydir <- getwd()   # tempdir()
## assume lead site1: enter "1" to allow transferring the control file
pda(site_id = 'site1', control = control, dir = mydir, upload_without_confirm = T,silent_message =F)

# S=readline(prompt="Type  <Return>   to continue : ")
## assume remote site3: enter "1" to allow tranferring your local estimate
pda(site_id = 'site3', ipdata = data_split[[3]], dir=mydir, upload_without_confirm = T)
  
# S=readline(prompt="Type  <Return>   to continue : ")
## assume remote site2: enter "1" to allow tranferring your local estimate
pda(site_id = 'site2', ipdata = data_split[[2]], dir=mydir, upload_without_confirm = T)

# S=readline(prompt="Type  <Return>   to continue : ")
## assume lead site1: enter "1" to allow tranferring your local estimate
## control.json is also automatically updated
pda(site_id = 'site1', ipdata = data_split[[1]], dir=mydir, upload_without_confirm = T)

# S=readline(prompt="Type  <Return>   to continue : ")
# ###################   STEP 2: derivative  #################
## assume remote site3: enter "1" to allow tranferring your derivatives
pda(site_id = 'site3', ipdata = data_split[[3]], dir=mydir,upload_without_confirm=T)

# S=readline(prompt="Type  <Return>   to continue : ")
## assume remote site2: enter "1" to allow tranferring your derivatives
pda(site_id = 'site2', ipdata = data_split[[2]], dir=mydir,upload_without_confirm=T)

# S=readline(prompt="Type  <Return>   to continue : ")
## assume lead site1: enter "1" to allow tranferring your derivatives
pda(site_id = 'site1', ipdata = data_split[[1]], dir=mydir,upload_without_confirm=T)

# S=readline(prompt="Type  <Return>   to continue : ")
# ############################  STEP 3: estimate  ###############################
## assume lead site1: enter "1" to allow tranferring the surrogate estimate
pda(site_id = 'site1', ipdata = data_split[[1]], dir=mydir,upload_without_confirm=T)

## all coef est:
# pooled:                     [-0.5129684,  0.1277404, 0.2421812, -0.2374705, 0.1359541
# meta:                       [-0.51927745,0.17385216,0.25506044,-0.25316367,-0.04065175]
# median:                     [-0.51845429, 0.23836445,0.61335488,-0.50068871,0]
# weighted.median:            [-0.51845429, 0.23836445,0.61335488,-0.50068871,0]
# lead:                       [-0.51845429, 0.48055952,0.89290808,-0.50068871,-0.17153903]
# odach_cc(meta):             [-0.51303977,0.12780133,0.24209696,-0.23769306,0.137274]
# odach_cc(median):           [-0.51369825, 0.11965385,0.24015613,-0.23640535,0.13420135]
# odach_cc(weighted.median):  [-0.51369825, 0.11965385,0.24015613,-0.23640535,0.13420135]
# odach_cc(lead):             [-0.5158358,  0.08234449,0.21134127,-0.24025793,0.13681052]

## all var est:
# "setilde": [0.11504608,0.35832839,0.36592225,0.31974574,0.34950701]
# cch_pooled  0.1150233, 0.3576904, 0.3655076, 0.3193014, 0.3499519 

# btilde and setilde (robust) are almost identical to pooled cch! 


# if you want to check ODACH_CC s.e. from inv-Hessian...
# config <- getCloudConfig(site_id =control$lead_site, dir=mydir)
# fit.pda = pdaGet(name = 'site1_estimate', config = config)
# sqrt(diag(solve(fit.pda$Htilde)))  


## the PDA ODACH_CC is now completed!
## All the sites can still run their own surrogate estimates and broadcast them.
 