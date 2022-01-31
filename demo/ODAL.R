require(survival)
require(data.table)
require(pda)
require(imager)
data(lung)

## In the toy example below we aim to analyze the association of lung status with age and sex using logistic regression,
## data(lung) from 'survival', we randomly assign to 3 sites: 'site1', 'site2', 'site3'
## we demonstrate using PDA ODAL can obtain a surrogate estimator that is close to the pooled estimate. 
## We run the example in local directory. In actual collaboration, account/password for pda server 
## will be assigned to the sites at the server https://pda.one.
## Each site can access via web browser to check the communication of the summary stats.

# Create 3 sites, split the lung data amongst them
sites = c('site1', 'site2', 'site3')
set.seed(42)
lung2 <- lung[,c('status', 'age', 'sex')]
lung2$sex <- lung2$sex - 1
lung2$status <- ifelse(lung2$status == 2, 1, 0)
lung_split <- split(lung2, sample(1:length(sites), nrow(lung), replace=TRUE))
## fit logistic reg using pooled data
fit.pool <- glm(status ~ age + sex, family = 'binomial', data = lung2)


S=readline(prompt="Type  <Return>   to continue : ")
# ############################  STEP 1: initialize  ###############################
control <- list(project_name = 'Lung cancer study',
                step = 'initialize',
                sites = sites,
                heterogeneity = FALSE,
                model = 'ODAL',
                family = 'binomial',
                outcome = "status",
                variables = c('age', 'sex'),
                # xlev = list(sex=c('F', 'M')),  #levels of all categorical X's, with the first being the reference
                optim_maxit = 100,
                lead_site = 'site1',
                upload_date = as.character(Sys.time()) )
 
 
## run the example in local directory:
## specify your working directory, default is the tempdir
mydir <- getwd()  # tempdir()
## assume lead site1: enter "1" to allow transferring the control file 
pda(site_id = 'site1', control = control, dir = mydir)
plot(load.image("https://raw.githubusercontent.com/Penncil/pda/master/demo/figures/p0_1.png"), axes=FALSE)

 
## in actual collaboration, account/password for pda server will be assigned, thus:
# pda(site_id = 'site1', control = control, uri = 'https://pda.one', secret='abc123')
## you can also set your environment variables, and no need to specify them in pda:
# Sys.setenv(PDA_USER = 'site1', PDA_SECRET = 'abc123', PDA_URI = 'https://pda.one')
# pda(site_id = 'site1', control = control)

S=readline(prompt="Type  <Return>   to continue : ")
## assume remote site3: enter "1" to allow tranferring your local estimate 
pda(site_id = 'site3', ipdata = lung_split[[3]], dir=mydir)
plot(load.image("https://raw.githubusercontent.com/Penncil/pda/master/demo/figures/p12_1.png"), axes=FALSE)

S=readline(prompt="Type  <Return>   to continue : ")
## assume remote site2: enter "1" to allow tranferring your local estimate  
pda(site_id = 'site2', ipdata = lung_split[[2]], dir=mydir)
plot(load.image("https://raw.githubusercontent.com/Penncil/pda/master/demo/figures/p11_1.png"), axes=FALSE)

S=readline(prompt="Type  <Return>   to continue : ")
## assume lead site1: enter "1" to allow tranferring your local estimate  
## control.json is also automatically updated
pda(site_id = 'site1', ipdata = lung_split[[1]], dir=mydir)
plot(load.image("https://raw.githubusercontent.com/Penncil/pda/master/demo/figures/p13_1.png"), axes=FALSE)

## if lead site1 initialized before other sites,
## lead site1: uncomment to sync the control before STEP 2
# pda(site_id = 'site1', control = control)
# config <- getCloudConfig(site_id = 'site1')
# pdaSync(config)

S=readline(prompt="Type  <Return>   to continue : ")
# ############################  STEP 2: derivative  ###############################
## assume remote site3: enter "1" to allow tranferring your derivatives  
pda(site_id = 'site3', ipdata = lung_split[[3]], dir=mydir)
plot(load.image("https://raw.githubusercontent.com/Penncil/pda/master/demo/figures/p22_1.png"), axes=FALSE)

S=readline(prompt="Type  <Return>   to continue : ")
## assume remote site2: enter "1" to allow tranferring your derivatives  
pda(site_id = 'site2', ipdata = lung_split[[2]], dir=mydir)
plot(load.image("https://raw.githubusercontent.com/Penncil/pda/master/demo/figures/p21_1.png"), axes=FALSE)

S=readline(prompt="Type  <Return>   to continue : ")
## assume lead site1: enter "1" to allow tranferring your derivatives  
pda(site_id = 'site1', ipdata = lung_split[[1]], dir=mydir)
 

S=readline(prompt="Type  <Return>   to continue : ")
# ############################  STEP 3: estimate  ###############################
## assume lead site1: enter "1" to allow tranferring the surrogate estimate  
pda(site_id = 'site1', ipdata = lung_split[[1]], dir=mydir)
plot(load.image("https://raw.githubusercontent.com/Penncil/pda/master/demo/figures/p31_1.png"), axes=FALSE)
 
## the PDA ODAL is now completed!
## All the sites can still run their own surrogate estimates and broadcast them.

S=readline(prompt="Type  <Return>   to continue : ")
## compare the surrogate estimate with the pooled estimate
config <- getCloudConfig(site_id = 'site1', dir=mydir)
fit.odal <- pdaGet(name = 'site1_estimate', config = config)
cbind(b.pool=fit.pool$coef,
      b.odal=fit.odal$btilde,
      sd.pool=summary(fit.pool)$coef[,2],
      sd.odal=sqrt(diag(solve(fit.odal$Htilde)/nrow(lung2))))

# S=readline(prompt="Type  <Return>   to continue : ")
# ## assume remote site2: (optional)
# pda(site_id = 'site2', ipdata = lung_split[[2]], dir=mydir)
# 
# S=readline(prompt="Type  <Return>   to continue : ")
# ## assume remote site3: (optional)
# pda(site_id = 'site3', ipdata = lung_split[[3]], dir=mydir)
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
# pda(site_id = 'site1', ipdata = lung_split[[1]], dir=mydir)
