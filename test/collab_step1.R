# rm(list=ls())
require(data.table)
source('R/PDA_engine.R')
mydata = fread(paste0('data/Lung_',Sys.getenv('PDA_SITE'),'.csv'))
## put intialize data on webdav
pda(data = mydata)
## waiting for other sites to initialize, and master site update pda_control
