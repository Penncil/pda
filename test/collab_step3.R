# rm(list=ls())
## (optional)  further synthesize 
require(data.table)
source('R/PDA_engine.R')
mydata = fread(paste0('data/Lung_',Sys.getenv('PDA_SITE'),'.csv'))
b_surr_syn = pda(data = mydata)
c(b_surr_syn$btilde)
