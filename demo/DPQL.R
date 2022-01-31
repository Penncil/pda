require(lme4)
require(minqa)
require(data.table)
require(pda) 
## In the toy example below we aim to rank the sites regarding their covid mortality
## after adjusting patients' age, sex and lab test using GLMM, 
## data(covid) is simulated and assumed to come from 6 sites,
## we demonstrate using PDA DPQL can obtain identical estimation as the pooled analysis.
## We run the example in local directory. In actual collaboration, account/password for pda server 
## will be assigned to the sites at the server https://pda.one.
## Each site can access via web browser to check the communication of the summary stats.

data(covid)  

# use site 1 2 3 as example
sites <- c('site1', 'site2', 'site3')
covid <- covid[site%in%sites,]
dd <- split(covid, covid$site)

S=readline(prompt="Type  <Return>   to continue : ")
# ############################ set DPQL control ###############################
control <- list(project_name = 'Covid mortality hospital profiling',
                step = 'derive_1',   # use zero as init thus no initialize step
                sites = sites,
                heterogeneity = T, 
                model = 'DPQL',
                family = 'binomial',
                outcome = "death",
                variables = c('age', 'sex', 'lab'),
                # xlev = list(sex=c('F', 'M')),  #levels of all categorical X's, with the first being the reference
                variables_heterogeneity = c('Intercept'),
                round = 1,
                maxround = 5,
                optim_maxit = 100,
                lead_site = 'site1',
                upload_date = as.character(Sys.time()) )


## run the example in local directory:
## specify your working directory, default is the tempdir
mydir <- getwd()  # tempdir()
## assume lead site1: enter "1" to allow transferring the control file 
pda(site_id = 'site1', control = control, dir = mydir)

# iterate deriv (collaborating sites calculate AD) and estimate (lead site) steps within each round
for(ir in 1:control$maxround){
  S=readline(prompt="Type  <Return>   to continue : ")
  ## assume remote site3: enter "1" to allow tranferring your aggregated data  
  pda(site_id = 'site3', ipdata = dd[[3]], dir=mydir)
  
  S=readline(prompt="Type  <Return>   to continue : ")
  ## assume remote site2: enter "1" to allow tranferring your AD
  pda(site_id = 'site2', ipdata = dd[[2]], dir=mydir)
  
  S=readline(prompt="Type  <Return>   to continue : ")
  ## assume lead site1: enter "1" to allow tranferring your Ad  
  ## control.json is also automatically updated
  pda(site_id = 'site1', ipdata = dd[[1]], dir=mydir)
  
  ## assume lead site1: enter "1" to allow tranferring the intermediate estimate  
  pda(site_id = 'site1', ipdata = dd[[1]], dir=mydir)
}
 

## the PDA DPQL is now completed! 
## compare the DPQL estimate with the pooled PQL estimate 
fit.dpql <- pdaGet(name = paste0('site1_estimate_',control$maxround), 
                   config = getCloudConfig(site_id = 'site1', dir=mydir))

fit.pool <- glmmPQL(death~age+sex+lab, ~1|site, data=covid, family='binomial') # 5 iteration
 
# fixef and ranef
cbind(bu.pool=c(fixef(fit.pool), ranef(fit.pool)$`(Intercept)`),
      bu.dpql=c(fit.dpql$bhat,fit.dpql$uhat) )
 
# var component
c(as.numeric(VarCorr(fit.pool)[1,1]), fit.dpql$Vhat) 
 
