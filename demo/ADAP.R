require(data.table)
require(pda)
require(glmnet)
data(ADAP_data) # load simulated dataset 

## In the toy example below we aim to analyze the association of a binary outcome 'status' 
## with 49 covariates 'x1' to 'x49'  using logistic lasso regression.
## In this toy example, we have 3 sites: 'site1', 'site2', 'site3' and each has 300 subjects.
## We demonstrate using PDA ADAP can obtain a surrogate estimator that is close to the pooled estimate. 
## We run the example in local directory. In actual collaboration, account/password for pda server 
## will be assigned to the sites at the server https://pda.one.
## Each site can access via web browser to check the communication of the summary stats.

ADAP_data <- data.frame(sites=ADAP_data$sites, status=ADAP_data$y, x=ADAP_data$x)
colnames(ADAP_data) <- c("sites", "status", paste0("x", 1:49))

## fit logistic reg using pooled data
fit.pool <- cv.glmnet(as.matrix(ADAP_data[,-c(1:2)]), ADAP_data$status, family = "binomial")
beta.pool <- as.matrix(coef(fit.pool, s = "lambda.min"))


sites <- unique(ADAP_data$sites)


S=readline(prompt="Type  <Return>   to continue : ")
# ############################  STEP 1: initialize  ###############################
control <- list(project_name = 'ADAP example',
                step = 'initialize',
                sites = sites,
                heterogeneity = FALSE,
                model = 'ADAP',
                family = 'lasso',
                outcome = "status",
                variables = colnames(ADAP_data)[-c(1,2)], 
                optim_maxit = 100, # not used
                lead_site = 'site1',
                upload_date = as.character(Sys.time()) )


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
pda(site_id = 'site3', ipdata = ADAP_data[sites == 'site3', -1], dir=mydir)

S=readline(prompt="Type  <Return>   to continue : ")
## assume remote site2: enter "1" to allow tranferring your local estimate  
pda(site_id = 'site2', ipdata = ADAP_data[sites == 'site2', -1], dir=mydir)

S=readline(prompt="Type  <Return>   to continue : ")
## assume lead site1: enter "1" to allow tranferring your local estimate  
## control.json is also automatically updated
pda(site_id = 'site1', ipdata = ADAP_data[sites == 'site1', -1], dir=mydir)

## if lead site1 initialized before other sites,
## lead site1: uncomment to sync the control before STEP 2
# pda(site_id = 'site1', control = control)
# config <- getCloudConfig(site_id = 'site1')
# pdaSync(config)

S=readline(prompt="Type  <Return>   to continue : ")
# ############################  STEP 2: derivative  ###############################
## assume remote site3: enter "1" to allow tranferring your derivatives  
pda(site_id = 'site3', ipdata = ADAP_data[sites == 'site3', -1], dir=mydir)

S=readline(prompt="Type  <Return>   to continue : ")
## assume remote site2: enter "1" to allow tranferring your derivatives  
pda(site_id = 'site2', ipdata = ADAP_data[sites == 'site2', -1], dir=mydir)

S=readline(prompt="Type  <Return>   to continue : ")
## assume lead site1: enter "1" to allow tranferring your derivatives  
pda(site_id = 'site1', ipdata = ADAP_data[sites == 'site1', -1], dir=mydir)
 
S=readline(prompt="Type  <Return>   to continue : ")
# ############################  STEP 3: estimate  ###############################
## assume lead site1: enter "1" to allow tranferring the surrogate estimate  
pda(site_id = 'site1', ipdata = ADAP_data[sites == 'site1', -1], dir=mydir)

## the PDA ADAP is now completed!

S=readline(prompt="Type  <Return>   to continue : ")
## compare the surrogate estimate with the pooled estimate
config <- getCloudConfig(site_id = 'site1', dir=mydir)
fit.adap <- pdaGet(name = 'site1_estimate', config = config)
cbind(b.pool=c(beta.pool), b.adap=fit.adap$btilde)
 
