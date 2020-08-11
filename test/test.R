require(survival)
require(data.table)
require(Rcpp)
# require(ODACO)
# rm(list=ls())


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
# write.csv(lung2, file='../Lung.csv', row.names = F)
 
write.csv(lung2[inst=='site1', -'inst'], file='data/site1/Lung_site1.csv', row.names = F)
write.csv(lung2[inst=='site2', -'inst'], file='data/site2/Lung_site2.csv', row.names = F)
write.csv(lung2[inst=='site3', -'inst'], file='data/site3/Lung_site3.csv', row.names = F)


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
                    cloud_website = 'webdav')

sink("data/pda_control.txt")
print(pda_control)
sink()
saveRDS(pda_control, 'data/pda_control.RDS')
################################## END: setup the ODAL test  ###################################################






# rm(list=ls())
mysite <- 'site1'    # provide your site abbre
mycloud <- "data/"   # provide your Dropbox PDA cloud folder 

setwd(mycloud)
list.files()

## eventually will be loading a package, 
## require(PDA)
## Rcpp::sourceCpp('../../PDA/src/rcpp_coxph.cpp')
source('../R/PDA_engine.R')

# pda_control = readRDS('pda_control.RDS')


setwd('../site1') 
mydata = fread('Lung_site1.csv')
pda(data = mydata, mysite='site1', mycloud=mycloud)
# -0.46925660    0.03174699   -1.52156716


## waiting for other sites to initialize


## master site update pda_control by adding beta_init (meta estimate), and step=2
pda_control_update(mycloud=mycloud)
# readRDS('pda_control.RDS')


## derivatives
setwd('../site1')
mydata = fread('Lung_site1.csv')
pda(data = mydata, mysite='site1', mycloud=mycloud)
 

## waiting for other sites to upload derivatives


## master site update pda_control$step = 3
pda_control_update(mycloud=mycloud)



## surrogate_est
setwd('../site1')
mydata = fread('Lung_site1.csv')
b_surr = pda(data = mydata, mysite='site1', mycloud=mycloud)


 
## compare
lung2 = fread('../Lung.csv')
fit.pool <- glm(status ~ age + sex, data = lung2, family = "binomial")
fit.pool$coef
rbind(pooled= fit.pool$coef,
      meta = pda_control$beta_init,
      odal = b_surr$btilde) 



## (optional) 
## master site update pda_control$step = 4 for further synthesize 
pda_control_update(mycloud=mycloud)

setwd('../site1')
# mydata = fread('Lung_site3.csv')
b_surr_syn = pda(data = mydata, mysite='site1', mycloud=mycloud)

c(b_surr_syn$btilde)







 
