require(survival)
require(data.table)
require(Rcpp)
################################## setup the test ODAL ###################################################
## prepare test data  
data(lung)
lung2=lung[,1:5]
lung2$inst[lung$inst>=2] = 2
lung2$inst[lung$inst>=12] = 3
lung2$inst[is.na(lung2$inst)] <- 3
lung2$inst <- c('site2', 'site3', 'site1')[lung2$inst]   # LETTERS[lung2$inst]
table(lung2$inst)
lung2$status <- ifelse(lung2$status == 2, 1, 0)
lung2$sex <- lung2$sex-1
lung2 <- data.table(lung2)
write.csv(lung2[inst=='site1', -'inst'], file='data/Lung_site1.csv', row.names = F)
write.csv(lung2[inst=='site2', -'inst'], file='data/Lung_site2.csv', row.names = F)
write.csv(lung2[inst=='site3', -'inst'], file='data/Lung_site3.csv', row.names = F)
write.csv(lung2[,-'inst'], file='data/Lung.csv')
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
mysite <- 'site1'    # provide your site abbre
## eventually will be loading a package, 
## require(PDA)
## Rcpp::sourceCpp('../../PDA/src/rcpp_coxph.cpp')
source('R/PDA_engine.R')
pda_put(pda_control,'pda_control',mysite)
# pda_control = readRDS('pda_control.RDS')
mydata = fread('data/Lung_site1.csv')
pda(data = mydata, mysite='site1')
# -0.46925660    0.03174699   -1.52156716
## waiting for other sites to initialize
## master site update pda_control by adding beta_init (meta estimate), and step=2
pda_control_update()
## derivatives
pda(data = mydata, mysite='site1')
## waiting for other sites to upload derivatives
## master site update pda_control$step = 3
pda_control_update()
## surrogate_est
b_surr = pda(data = mydata, mysite='site1')
## compare
lung2 = fread('data/Lung.csv')
fit.pool <- glm(status ~ age + sex, data = lung2, family = "binomial")
fit.pool$coef
rbind(pooled= fit.pool$coef,
      meta = pda_control$beta_init,
      odal = b_surr$btilde) 
## (optional) 
## master site update pda_control$step = 4 for further synthesize 
pda_control_update()
b_surr_syn = pda(data = mydata, mysite='site1')
c(b_surr_syn$btilde)
