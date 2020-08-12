# rm(list=ls())
require(data.table)
source('R/PDA_engine.R')
mydata = fread(paste0('data/Lung_',Sys.getenv('PDA_SITE'),'.csv'))
#get pda_control data from server
pda_control<-pda_get('pda_control')
## surrogate_est
b_surr = pda(data = mydata)
## compare
lung2 = fread('data/Lung.csv')
fit.pool <- glm(status ~ age + sex, data = lung2, family = "binomial")
fit.pool$coef
rbind(pooled= fit.pool$coef,
      meta = pda_control$beta_init,
      odal = b_surr$btilde) 
