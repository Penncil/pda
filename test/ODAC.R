#
require(survival)
require(data.table)
require(pda)
data(lung)
source('R/pda.R')
source('R/ODAL.R')
#create 3 sites and split the lung data amongst them
all_site = c('site1', 'site2', 'site3')
set.seed(42)
lung2<-lung[,2:5]
lung2$sex <- lung2$sex-1
lung2$status <- ifelse(lung2$status == 2, 1, 0)
lung2<-split(lung2, sample(1:length(all_site), nrow(lung), replace=T))
############################### setup  ODAL ##############################
## RUN BY LEAD ONLY 
control <- list(project_name = 'Lung cancer study',
        step = 'initialize',      # current step of iteration. updated by lead site  
        all_site = all_site,
        common_data_model = 'OMOP',
        heterogeneity = FALSE,
        model = 'ODAL',
        family = 'binomial',
        outcome = "status",
        variables = c('age','sex'),
        lead_site = all_site[1],
        last_update = as.character(Sys.time()),
        upload_date = as.Date('2020-08-15'),
        optim_maxit=100,
        reference = 'https://academic.oup.com/jamia/article-abstract/27/3/376/5670808',
        heterogeneity = FALSE)
pda_create(control)
################################## END: setup the ODAL test  ################
#run pda to seed data on server
#initialize
while (is.character(control$step)) {
  print(paste("step:",step))
  for(i in 1:length(all_site)) {
    pda(ipdata=lung2[[i]],site_id=all_site[i])
  }
  control<-control_update()
}
