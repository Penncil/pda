require(lme4)
require(minqa)
require(data.table)
require(pda) 
## In the toy example below we aim to analyze the association of hospitalization length of stay (LOS) 
## with age, sex and lab test using linear model, 
## data(LOS) is simulated and assumed to come from 3 sites: 'site1', 'site2', 'site3'
## we demonstrate using PDA DLM can obtain identical estimation as the pooled analysis.
## 3 models are tested:
# - Linear model ignoring site-effect
# - Linear model with fixed site-effect
# - Linear model with random site-effect (Linear mixed model)
## We run the example in local directory. In actual collaboration, account/password for pda server 
## will be assigned to the sites at the server https://pda.one.
## Each site can access via web browser to check the communication of the summary stats.

data(LOS)  
sites = c('site1', 'site2', 'site3')
LOS_split <- split(LOS, LOS$site)

S=readline(prompt="Type  <Return>   to continue : ")
# ############################  STEP 1: initialize  ###############################
control <- list(project_name = 'Length of stay study',
                step = 'initialize',
                sites = sites,
                heterogeneity = FALSE,
                # heterogeneity_effect = 'fixed',
                model = 'DLM',
                family = 'gaussian',
                outcome = "los",
                variables = c('age', 'sex', 'lab'),
                # xlev = list(sex=c('F', 'M')),  #levels of all categorical X's, with the first being the reference
                # variables_heterogeneity = c('Intercept'),
                optim_maxit = 100,
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
pda(site_id = 'site3', ipdata = LOS_split[[3]], dir=mydir)
 
S=readline(prompt="Type  <Return>   to continue : ")
## assume remote site2: enter "1" to allow tranferring your local estimate  
pda(site_id = 'site2', ipdata = LOS_split[[2]], dir=mydir)
 
S=readline(prompt="Type  <Return>   to continue : ")
## assume lead site1: enter "1" to allow tranferring your local estimate  
## control.json is also automatically updated
pda(site_id = 'site1', ipdata = LOS_split[[1]], dir=mydir)
 

S=readline(prompt="Type  <Return>   to continue : ")
# ############################  STEP 2: estimate  ###############################
## assume lead site1: enter "1" to allow tranferring the surrogate estimate  
pda(site_id = 'site1', ipdata = LOS_split[[1]], dir=mydir)
 
## the PDA DLM is now completed!

S=readline(prompt="Type  <Return>   to continue : ")
## compare the DLM estimate with the pooled estimate 
fit.pool <- lm(los~age+sex+lab, data=LOS)  # LM ignoring site effect

config <- getCloudConfig(site_id = 'site1', dir=mydir)
fit.dlm <- pdaGet(name = 'site1_estimate', config = config)
cbind(b.pool=fit.pool$coef,
      b.dlm=fit.dlm$bhat,
      sd.pool=summary(fit.pool)$coef[,2],
      sd.dlm=fit.dlm$sebhat)

 
S=readline(prompt="Type  <Return>   to continue : ")
# # ########################  DLM with fixed site-effect  ###################### 
# reset control to re-fit DLM with fixed site-effect, no extra AD communication is needed 
# Jessie, please double check the below demo to make sure it works
control <- list(project_name = 'Length of stay study',
                step = 'estimate',
                sites = sites,
                heterogeneity = TRUE,
                heterogeneity_effect = 'fixed',
                model = 'DLM',
                family = 'gaussian',
                outcome = "los",
                variables = c('age', 'sex', 'lab'),
                # variables_heterogeneity = c('Intercept'),
                optim_maxit = 100,
                lead_site = 'site1',
                upload_date = as.character(Sys.time()) )
pda(site_id = 'site1', control = control, dir = mydir)
pda(site_id = 'site1', ipdata = LOS_split[[1]], dir=mydir)

## fit LM using pooled data, assuming fixed site effect
fit.pool <- lm(los~age+sex+lab+site, data=LOS)
# config <- getCloudConfig(site_id = 'site1', dir=mydir)
fit.dlm <- pdaGet(name = 'site1_estimate', config = config)
cbind(b.pool=fit.pool$coef,
      b.dlm=c(fit.dlm$bhat, fit.dlm$uhat),       
      sd.pool=summary(fit.pool)$coef[,2],  
      sd.dlm=c(fit.dlm$sebhat, fit.dlm$seuhat))  




# # ########################  DLM with random site-effect (DLMM) ###################### 
# reset control to re-fit DLMM, no extra AD communication is needed 
control <- list(project_name = 'Length of stay study',
                step = 'estimate',
                sites = sites,
                heterogeneity = TRUE,
                heterogeneity_effect = 'random',
                model = 'DLM',
                family = 'gaussian',
                outcome = "los",
                variables = c('age', 'sex', 'lab'),
                # variables_heterogeneity = c('Intercept'),
                optim_maxit = 100,
                lead_site = 'site1',
                upload_date = as.character(Sys.time()) )
pda(site_id = 'site1', control = control, dir = mydir)
pda(site_id = 'site1', ipdata = LOS_split[[1]], dir=mydir)

## fit LM using pooled data, assuming fixed site effect
fit.pool <- lme4::lmer(los~age+sex+lab+(1|site), REML = F, data=LOS)
fit.dlm <- pdaGet(name = 'site1_estimate', config = config)
# fixed effects
cbind(b.pool = summary(fit.pool)$coef[,1],
      b.dlm = c(fit.dlm$bhat),      
      sd.pool = summary(fit.pool)$coef[,2],  
      sd.dlm = fit.dlm$sebhat)  

# var component
cbind(data.frame(summary(fit.pool)$varcor)$vcov, 
      c(fit.dlm$Vhat, fit.dlm$sigmahat^2) )

# random effects (BLUP)
cbind(u.pool = ranef(fit.pool)$site,
      u.dlm = c(fit.dlm$uhat))
