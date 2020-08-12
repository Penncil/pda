require(data.table)
require(Rcpp)
source('R/PDA_engine.R')
## master site update pda_control by adding beta_init (meta estimate), and step=2
pda_control_update()
## derivatives
mydata = fread('data/Lung_site1.csv')
pda(data = mydata)
## waiting for other sites to upload derivatives
