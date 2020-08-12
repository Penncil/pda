# rm(list=ls())
## (optional)  further synthesize 
require(data.table)
source('R/PDA_engine.R')
sitename<-Sys.getenv('PDA_SITE')
mydata = fread(paste0('data/Lung_',sitename,'.csv'))
b_surr_syn = pda(data = mydata)
c(b_surr_syn$btilde)
