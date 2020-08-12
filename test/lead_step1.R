require(data.table)
require(Rcpp)
source('R/PDA_engine.R')
################################## setup the test ODAL ###################################################
## research project setting, decided by all collaborators before analysis 
pda_control <- list(project_name = 'Lung cancer study',
                    step = 1,      # current step of iteration. updated by the master site  
                    all_site = c('site1', 'site3', 'site2'),
                    common_data_model = 'OMOP',
                    formula = status ~ age + sex,
                    family = 'binomial',
                    heterogeneity = FALSE,
                    PDA_model = 'ODAL', 
                    master_site = 'site1',
                    last_update = as.character(Sys.time()),
                    upload_date = as.Date('2020-08-15'),
                    optim_maxit=100,
                    reference = 'https://academic.oup.com/jamia/article-abstract/27/3/376/5670808',
                    heterogeneity = FALSE)
################################## END: setup the ODAL test  ###################################################
pda_put(pda_control,'pda_control')
sitename<-Sys.getenv('PDA_SITE')
mydata = fread(paste0('data/Lung_',sitename,'.csv'))
#run pda to seed data on server
pda(data = mydata)
# -0.46925660    0.03174699   -1.52156716
## waiting for other sites to initialize

