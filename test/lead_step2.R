require(data.table)
require(Rcpp)
source('R/PDA_engine.R')
sitename<-Sys.getenv('PDA_SITE')
mydata = fread(paste0('data/Lung_',sitename,'.csv'))
## master site update pda_control by adding beta_init (meta estimate), and step=2
pda_control_update()
## derivatives
pda(data = mydata)
## waiting for other sites to upload derivatives
