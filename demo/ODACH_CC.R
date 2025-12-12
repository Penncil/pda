require(survival)
require(data.table)
# require(pda) 

## In the toy example below we aim to analyze the association of survival {time} 
##  with X1, Category, and Group using Cox reg with case-cohort design,
## data(odach_cc) simulated, subsampled case+cohort data, 3 sites: 'site1', 'site2', 'site3'
## we demonstrate using PDA ODACH_CC can obtain a surrogate estimator that is close to the pooled estimate.
## ODACH_CC relies on the surrogate of the pooled stratified case-cohort pseudo likelihood function
## We run the example in local directory. 
## In actual collaboration, the data communication can be done via the PDA_OTA platform https://pda-ota.pdamethods.org/
## Each site can access via web browser to transfer aggregate data and check the progress of the project.


# load('/Users/chongliang/Dropbox/PDA-git/pda/data/odach_cc.rda')
data(odach_cc)
odach_cc$Category = as.character(odach_cc$Category)
# to test site3 has missing category Z for "Category" variable
odach_cc[odach_cc$site=='site3'&odach_cc$Category=='Z (X,Y,Z)', ]$Category='Y (X,Y,Z)'
table(odach_cc$site, odach_cc$Category)

# add time_in for time-interval outcome
precision <- min(diff(sort(odach_cc$time)))/2 # 1e-4
odach_cc$time_in = 0 # start time
odach_cc[odach_cc$subcohort == 0, "time_in"] <- odach_cc[odach_cc$subcohort == 0, "time"] - precision
odach_cc$ID = 1:nrow(odach_cc)
data_split <- split(odach_cc, odach_cc$site)


# pooled estimates
fit.pool1 <- coxph(Surv(time_in, time, status) ~ X1+X2 +Category+ strata(site,Group)+cluster(ID), data=odach_cc, robust=T)
fit.pool1$coef
#         X1                X2 CategoryY (X,Y,Z) CategoryZ (X,Y,Z) 
# -0.5412024         0.3930600        -0.1212722         0.1421058  
summary(fit.pool1)$coef[,'se(coef)']
# 0.1101166         0.1334716         0.2890987         0.3253451 


# ODACH_CC (DRAFT)
# ############################  STEP 1: initialize  ############################### 
sites = c('site1', 'site2', 'site3')
control <- list(project_name = 'ODACH case-cohort toy example',
                step = 'initialize',
                sites = sites,
                heterogeneity = TRUE,
                model = 'ODACH_CC',
                family = 'cox',
                # time-interval as outcome, "Surv(time, status)" is also OK
                outcome = "Surv(time_in, time, status)", 
                variables = c('X1','X2', 'Category'),
                #levels of all categorical X's, with the first being the reference
                variables_lev = list(Category=c('X (X,Y,Z)','Y (X,Y,Z)','Z (X,Y,Z)')),
                strata_names = c('Group'), # strata within sites
                method = 'Prentice',                      
                init_method = 'meta', 
                optim_method = 'BFGS',
                optim_maxit = 100,
                lead_site = 'site1',
                upload_date = as.character(Sys.time()) )
 
## run the example in local directory:
## specify your working directory, default is the tempdir
mydir <- getwd()   # tempdir()
file.remove(list.files(mydir,full.names = T)[grepl('.json', list.files(mydir))])

pda(site_id = 'site1', control = control, dir = mydir, upload_without_confirm = T,silent_message =F) # 

pda(site_id = 'site3', ipdata = data_split[[3]], dir=mydir, upload_without_confirm = T,silent_message =F)
  
pda(site_id = 'site2', ipdata = data_split[[2]], dir=mydir, upload_without_confirm = T,silent_message =F)

# control.json is also automatically updated
pda(site_id = 'site1', ipdata = data_split[[1]], dir=mydir, upload_without_confirm = T,silent_message =F)

# ###################   STEP 2: derivative  #################
pda(site_id = 'site3', ipdata = data_split[[3]], dir=mydir,upload_without_confirm=T)

pda(site_id = 'site2', ipdata = data_split[[2]], dir=mydir,upload_without_confirm=T)

pda(site_id = 'site1', ipdata = data_split[[1]], dir=mydir,upload_without_confirm=T)

# ############################  STEP 3: estimate  ###############################
pda(site_id = 'site1', ipdata = data_split[[1]], dir=mydir,upload_without_confirm=T)
# "btilde":[-0.5412,0.393,-0.1215,0.1424]
# "beta_init":[-0.5561,0.3879,-0.1403,-0.0218]
# b_pool      -0.5412,0.3931,-0.1213,0.1421 

# "setilde":[0.1142,0.1576,0.3056,0.3865]
# se_pool    0.1101,0.1335,0.2891,0.3253 

 
 
# if you want to check ODACH_CC s.e. from inv-Hessian...
# config <- getCloudConfig(site_id =control$lead_site, dir=mydir)
# fit.pda = pdaGet(name = 'site1_estimate', config = config)
# sqrt(diag(solve(fit.pda$Htilde)))  


## the PDA ODACH_CC is now completed!