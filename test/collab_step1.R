# rm(list=ls())
require(data.table)
source('R/PDA_engine.R')
## initialize
sitename<-Sys.getenv('PDA_SITE')
mydata = fread(paste0('data/Lung_',sitename,'.csv'))
pda(data = mydata)
## waiting for other sites to initialize, and master site update pda_control
