require(survival)
require(data.table)
require(pda)
data(lung)

## Create a number of sites, split the lung data amongst them
sites = c('site1', 'site2', 'site3')
set.seed(42)
lung2 <- lung[,2:5]
lung2$sex <- lung2$sex - 1
lung2$status <- ifelse(lung2$status == 2, 1, 0)
lung_split <- split(lung2, sample(1:length(sites), nrow(lung), replace=TRUE))
## fit logistic reg using pooled data
fit.pool <- glm(status ~ age + sex, family = 'binomial', data = lung2)

## In the example below we aim to use PDA ODAL to obtain a surrogate estimator that is close to the pooled estimate
## Accounts (site1, site2, site3) and passwords (WLjySoaZnSqMNswowg) are given to the 3 example sites 
## at the server https://pda.one. Each site can access via web browser to check the communication of the summary stats.

# ############################  STEP 1: initialize  ###############################
## lead site1, please review and enter "1" to allow putting the project control file to the server
control <- list(project_name = 'Lung cancer study',
                step = 'initialize',     
                sites = sites,
                heterogeneity = FALSE,
                model = 'ODAL',
                outcome = "status",
                variables = c('age', 'sex'),
                optim_maxit = 100,
                lead_site = sites[1],
                upload_date = as.character(Sys.time()) )
Sys.setenv(PDA_USER = 'site1', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
pda(site_id = 'site1', control = control)  
# now the group would see control.json	at https://pda.one

# remote site3, please review and enter "1" to allow putting your local estimation to the server  
i <- 3
Sys.setenv(PDA_USER = 'site3', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
control <- pda(ipdata = lung_split[[i]], site_id = sites[i])
# now the group would see site3_initialize.json	at https://pda.one

## remote site2, please review and enter "1" to allow putting your local estimation to the server
i <- 2
Sys.setenv(PDA_USER = 'site2', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
control <- pda(ipdata = lung_split[[i]], site_id = sites[i])
# now the group would see site2_initialize.json	at https://pda.one

## lead site1, please review and enter "1" to allow putting your local estimation to the server
i <- 1
Sys.setenv(PDA_USER = 'site1', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
control <- pda(ipdata = lung_split[[i]], site_id = sites[i]) 
# now the group would see site1_initialize.json	at https://pda.one
# control.json is also automatically updated  

## if lead site1 initialized before other sites, 
## lead site1: uncomment to synchoronize the control before STEP 2
# Sys.setenv(PDA_USER = 'site1', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
# pda(site_id = 'site1', control = control)
# config <- getCloudConfig(site_id = 'site1')
# pdaSync(config)

# ############################  STEP 2: derivative  ###############################
## remote site3: please review and enter "1" to allow putting your derivatives to the server   
i <- 3
Sys.setenv(PDA_USER = 'site3', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
control <- pda(ipdata = lung_split[[i]], site_id = sites[i])
# now the group would see site3_derive.json	at https://pda.one

## remote site2: please review and enter "1" to allow putting your derivatives to the server    
i <- 2
Sys.setenv(PDA_USER = 'site2', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
control <- pda(ipdata = lung_split[[i]], site_id = sites[i])
# now the group would see site2_derive.json	at https://pda.one

## lead site1: please review and enter "1" to allow putting your derivatives to the server      
i <- 1
Sys.setenv(PDA_USER = 'site1', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
control <- pda(ipdata = lung_split[[i]], site_id = sites[i])
# now the group would see site1_derive.json	at https://pda.one

# ############################  STEP 3: estimate  ###############################
## lead site1: please review and enter "1" to allow putting the surrogate estimate to the server     
i <- 1
Sys.setenv(PDA_USER = 'site1', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
control <- pda(ipdata = lung_split[[i]], site_id = sites[i])
# now the group would see site1_estimate.json	at https://pda.one
## compare the surrogate estimate with the pooled estimate
fit.odal <- pdaGet(name = 'site1_estimate', config = config)
cbind(b.pool=fit.pool$coef, 
      b.odal=fit.odal$btilde, 
      sd.pool=summary(fit.pool)$coef[,2], 
      sd.odal=sqrt(diag(solve(fit.odal$Htilde)/nrow(lung2))))
## the PDA ODAL is now completed! All the sites can still run their surrogate estimates and broadcast them 

## remote site2: (optional)  
i <- 2
Sys.setenv(PDA_USER = 'site2', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
control <- pda(ipdata = lung_split[[i]], site_id = sites[i])

## remote site3: (optional)
i <- 3
Sys.setenv(PDA_USER = 'site3', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
control <- pda(ipdata = lung_split[[i]], site_id = sites[i]) 


## If all the sites broadcast their surrogate estimates, a final synthesize step will further improve the estimate.
## lead site1: uncomment to synchoronize the control before STEP 4
Sys.setenv(PDA_USER = 'site1', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
pda(site_id = 'site1', control = control)
config <- getCloudConfig(site_id = 'site1')
pdaSync(config)

# ########################  STEP 4: synthesize (optional)  ######################## 
## lead site1:     
i <- 1
Sys.setenv(PDA_USER = 'site1', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
control <- pda(ipdata = lung_split[[i]], site_id = sites[i])







Sys.setenv(PDA_USER = 'site1', 
           PDA_SECRET = 'WLjySoaZnSqMNswowg', 
           PDA_URI = 'https://pda.one/003')
pda(site_id = 'site1', control = control)

config <- getCloudConfig(site_id = 'site1')

# Run pda until step is empty
while (is.character(control$step)) {
  print(paste("step:", control$step))
  # Cycle through sites
  for(i in length(sites):1) {      # , PDA_DIR='/Users/chl18019/Dropbox/PDA-git'
    Sys.setenv(PDA_USER = paste0('site',i), 
               PDA_SECRET = 'WLjySoaZnSqMNswowg', 
               PDA_URI = 'https://pda.one/003')
    control<-pda(ipdata=lung_split[[i]],site_id=sites[i])
  }
}

# RUN BY LEAD ONLY , check results at https://pda.one/003 