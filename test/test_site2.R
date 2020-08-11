# rm(list=ls())
require(data.table)
mysite <- 'site2'     # provide your site abbre
# mycloud <- "/Users/chl18019/Dropbox/_PDAmethods/2_Software/Lung_ODAL/cloud/"   # provide your Dropbox PDA cloud folder 
source('R/PDA_engine.R')
# pda_control = readRDS('pda_control.RDS')
## initialize
mydata = fread('data/Lung_site2.csv')
pda(data = mydata, mysite='site2')
# -1.88570040    0.05577115   -1.06961062  
## waiting for other sites to initialize, and master site update pda_control
## derivatives
mydata = fread('data/Lung_site2.csv')
pda(data = mydata, mysite='site2')
## waiting for master site update pda_control$step = 3
## surrogate_est
mydata = fread('data/Lung_site2.csv')
b_surr = pda(data = mydata, mysite='site2')
## compare
lung2 = fread('data/Lung.csv')
fit.pool <- glm(status ~ age + sex, data = lung2, family = "binomial")
fit.pool$coef
rbind(pooled= fit.pool$coef,
      meta = pda_control$beta_init,
      odal = b_surr$btilde) 
## (optional)  further synthesize 
# mydata = fread('Lung_site2.csv')
b_surr_syn = pda(data = mydata, mysite='site2', mycloud=mycloud)
c(b_surr_syn$btilde)
