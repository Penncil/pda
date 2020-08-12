require(data.table)
require(Rcpp)
source('R/PDA_engine.R')
mydata = fread(paste0('data/Lung_',Sys.getenv('PDA_SITE'),'.csv'))
## master site update pda_control$step = 4 for further synthesize 
pda_control_update()
b_surr_syn = pda(data = mydata)
c(b_surr_syn$btilde)
