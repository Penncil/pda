require(survival)
require(data.table)
# require(pda) 
 
# load('/Users/chongliangluo/Library/CloudStorage/Dropbox/PDA-git/pda/data/odach_cc.rda')
# setwd('/Users/chongliangluo/Library/CloudStorage/Dropbox/PDA_test/CL/')

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
data_split <- split(odach_cc, odach_cc$site)

## cch stratified by site, using self-written cch_pooled()
fit.pool <- cch_pooled(Surv(time, status) ~ X1+`Group (A,B,C)`+Category, data=odach_cc,  
                       full_cohort_size=c(site1=800,site2=600,site3=400), method='Prentice')
fit.pool$par
#          X1   `Group (A,B,C)`B  `Group (A,B,C)`C CategoryY (X,Y,Z) CategoryZ (X,Y,Z) 
# -0.50824965         0.12309366        0.25814696      -0.25125968         0.03014395 
sqrt(diag(-solve(fit.pool$hessian)))
# 0.1097648         0.3183301         0.2966302         0.3044304         0.2829128 


# cch each site
sapply(1:3, function(i) cch(Surv(time, status) ~ X1+`Group (A,B,C)`+Category, 
                            data = cbind(ID=1:nrow(data_split[[i]]), data_split[[i]]),
                            subcoh = ~subcohort, id = ~ID, 
                            cohort.size =c(site1=800,site2=600,site3=400)[i], 
                            method = 'Prentice')$coef) 
#                          [,1]       [,2]         [,3]
# XX1                -0.5184543 -0.4707972 -0.610918824
# X`Group (A,B,C)`B   0.4805595  0.2383644 -0.074294738
# X`Group (A,B,C)`C   0.8929081  0.6133549 -1.371104126
# XCategoryY (X,Y,Z) -0.5006887 -0.7252045  0.747914104
# XCategoryZ (X,Y,Z) -0.1715390  0.1941511 -0.006210128

sites = c('site1', 'site2', 'site3')
S=readline(prompt="Type  <Return>   to continue : ")
# ############################  STEP 1: initialize  ###############################
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


## run the example in local directory:
## specify your working directory, default is the tempdir
mydir <- getwd()   # tempdir()
## assume lead site1: enter "1" to allow transferring the control file
pda(site_id = 'site1', control = control, dir = mydir, upload_without_confirm = T,silent_message =F)
## in actual collaboration, account/password for pda server will be assigned, thus:
# pda(site_id = 'site1', control = control, uri = 'https://pda.one', secret='abc123')
## you can also set your environment variables, and no need to specify them in pda:
# Sys.setenv(PDA_USER = 'site1', PDA_SECRET = 'abc123', PDA_URI = 'https://pda.one')
# pda(site_id = 'site1', control = control)

S=readline(prompt="Type  <Return>   to continue : ")
## assume remote site3: enter "1" to allow tranferring your local estimate
pda(site_id = 'site3', ipdata = data_split[[3]], dir=mydir, upload_without_confirm = T)
  
S=readline(prompt="Type  <Return>   to continue : ")
## assume remote site2: enter "1" to allow tranferring your local estimate
pda(site_id = 'site2', ipdata = data_split[[2]], dir=mydir, upload_without_confirm = T)

S=readline(prompt="Type  <Return>   to continue : ")
## assume lead site1: enter "1" to allow tranferring your local estimate
## control.json is also automatically updated
pda(site_id = 'site1', ipdata = data_split[[1]], dir=mydir, upload_without_confirm = T)


S=readline(prompt="Type  <Return>   to continue : ")
# ###################   STEP 2: derivative  #################
## assume remote site3: enter "1" to allow tranferring your derivatives
pda(site_id = 'site3', ipdata = data_split[[3]], dir=mydir,upload_without_confirm=T)

S=readline(prompt="Type  <Return>   to continue : ")
## assume remote site2: enter "1" to allow tranferring your derivatives
pda(site_id = 'site2', ipdata = data_split[[2]], dir=mydir,upload_without_confirm=T)

S=readline(prompt="Type  <Return>   to continue : ")
## assume lead site1: enter "1" to allow tranferring your derivatives
pda(site_id = 'site1', ipdata = data_split[[1]], dir=mydir,upload_without_confirm=T)

S=readline(prompt="Type  <Return>   to continue : ")
# ############################  STEP 3: estimate  ###############################
## assume lead site1: enter "1" to allow tranferring the surrogate estimate
pda(site_id = 'site1', ipdata = data_split[[1]], dir=mydir,upload_without_confirm=T)
# "btilde":[-0.50811268,0.12339569,0.25802334,-0.25067117,0.03074594]
# "setilde":[0.13333819,0.35353426,0.31805217,0.3313836,0.33198831]
# btilde is almost identical to pooled cch!
# setilde (robust) is slightly greater than inv-Hessian 

config <- getCloudConfig(site_id =control$lead_site, dir=mydir)
fit.pda = pdaGet(name = 'site1_estimate', config = config)
sqrt(diag(solve(fit.pda$Htilde))) # inv-Hessian s.e. est
# 0.1101771 0.3178109 0.2988028 0.3016307 0.2845887


## the PDA ODACH_CC is now completed!
## All the sites can still run their own surrogate estimates and broadcast them.
 