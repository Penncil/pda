# rm(list=ls())
require(data.table)

mysite <- 'site3'     # provide your site abbre
source('R/PDA_engine.R')


## initialize
mydata = fread('data/Lung_site3.csv')
pda(data = mydata, mysite='site3')
# 0.728449481   0.009317684  -0.296427795 


## waiting for other sites to initialize, and master site update pda_control


## derivatives
mydata = fread('data/site3/Lung_site3.csv')
pda(data = mydata, mysite='site3', mycloud=mycloud)


## waiting for master site update pda_control$step = 3


## surrogate_est
mydata = fread('data/site3/Lung_site3.csv')
b_surr = pda(data = mydata, mysite='site3', mycloud=mycloud)


## compare
lung2 = fread('data/Lung.csv')
fit.pool <- glm(status ~ age + sex, data = lung2, family = "binomial")
fit.pool$coef

rbind(pooled= fit.pool$coef,
      meta = pda_control$beta_init,
      odal = b_surr$btilde)


## (optional)  further synthesize 
# mydata = fread('Lung_site3.csv')
b_surr_syn = pda(data = mydata, mysite='site3', mycloud=mycloud)

c(b_surr_syn$btilde)
