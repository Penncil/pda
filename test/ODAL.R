require(survival)
require(data.table)
require(pda)
data(lung)
#create a number of sites, split the lung data amongst them
sites = c('site1', 'site2', 'site3')
set.seed(42)
lung2<-lung[,2:5]
lung2$sex <- lung2$sex-1
lung2$status <- ifelse(lung2$status == 2, 1, 0)
lung_split<-split(lung2, sample(1:length(sites), nrow(lung), replace=T))
######################### setup  ODAL #############################
control <- list(project_name = 'Lung cancer study',
        step = 'initialize',    # current step, updated by lead
        sites = sites,
        heterogeneity = FALSE,
        model = 'ODAL',
        outcome = "status",
        variables = c('age', 'sex'),
        optim_maxit=100,
        lead_site = sites[1],
        upload_date = as.character(Sys.time()),
        heterogeneity = FALSE)
## RUN BY LEAD ONLY 
pdaInit(control)
#run pda until step is empty
while (is.character(control$step)) {
  print(paste("step:",control$step))
  #cycle through sites
  for(i in 1:length(sites)) {
    pdaStep(ipdata=lung_split[[i]],site_id=sites[i])
  }
  #update control
  control<-pdaSync()
}
