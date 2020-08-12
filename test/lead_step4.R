require(data.table)
require(Rcpp)
## master site update pda_control$step = 4 for further synthesize 
pda_control_update()
b_surr_syn = pda(data = mydata)
c(b_surr_syn$btilde)
