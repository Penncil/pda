PDA: Privacy-preserving Distributed Algorithms
==============================================


## Outline

1. PDA workflow
2. Package requirements
3. Instructions for installing and running pda package
4. FAQ


## PDA Workflow
![](Picture1.png)

## Package Requirements
- A database with clear and consistent variable names
- On Windows: download and install [RTools](http://cran.r-project.org/bin/windows/Rtools/) 


## Instructions for Installing and Running pda Package

Below are the instructions for installing and then running the package.

### How to install the pda package
There are several ways in which one could install the `pda` package. 

1. [Q: do we need to setup the environment for C++ first?]

2. In RStudio, create a new project: File -> New Project... -> New Directory -> New Project. 

3. Execute the following R code: 

```r
# Install the latest version of PDA in R:
install.packages("pda")
library(pda)

# Or you can install via github:
install.packages("devtools")
library(devtools)
devtools::install_github("penncil/pda")
library(pda)
```

### How to run pda package

Below are two ways to run the pda example. ## In the example below we aim to use PDA ODAL to obtain a surrogate estimator that is close to the pooled estimate. 

#### *Run example with demo(pda)*

```r
demo(pda)
``` 

#### *Run example with code*

Step 0: load related R packages and prepare sample data

```r
# load packages
require(survival)
require(data.table)

# sample data, lung, from "survival" package
data(lung)

# create a number of sites, split the lung data amongst them
sites = c('site1', 'site2', 'site3')
set.seed(42)
lung2<-lung[,2:5]
lung2$sex <- lung2$sex-1
lung2$status <- ifelse(lung2$status == 2, 1, 0)
lung_split<-split(lung2, sample(1:length(sites), nrow(lung), replace=TRUE))
``` 
Step 1: Initialization

```r
## fit logistic reg using pooled data
fit.pool <- glm(status ~ age + sex, family = 'binomial', data = lung2)


# ############################  STEP 1: initialize  ###############################
## lead site1: please review and enter "1" to allow putting the control file to the server
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
# A cloud server will be used in a real collaborative project. The user name (i.e., PDA_USER) and the password (i.e., PDA_SECRET) will be assigned to the working sites by the lead site. The site with assigned username and code can upload the results to the cloud server. 
# Sys.setenv(PDA_USER = 'username', PDA_SECRET = 'password', PDA_URI = 'https://pda.one')
pda(site_id = 'site1', control = control)
# now the group would see control.json	at https://pda.one

## remote site2: please review and enter "1" to allow putting your local estimate to the server
i <- 2
# Sys.setenv(PDA_USER = 'site2', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
control <- pda(ipdata = lung_split[[i]], site_id = sites[i])
# now the group would see site2_initialize.json	at https://pda.one

# remote site3: please review and enter "1" to allow putting your local estimate to the server
i <- 3
# Sys.setenv(PDA_USER = 'site3', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
control <- pda(ipdata = lung_split[[i]], site_id = sites[i])
# now the group would see site3_initialize.json	at https://pda.one


## lead site1: please review and enter "1" to allow putting your local estimate to the server
i <- 1
# Sys.setenv(PDA_USER = 'site1', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
control <- pda(ipdata = lung_split[[i]], site_id = sites[i])
# now the group would see site1_initialize.json	at https://pda.one
# control.json is also automatically updated

## if lead site1 initialized before other sites,
## lead site1: uncomment to synchoronize the control before STEP 2
# Sys.setenv(PDA_USER = 'site1', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
# pda(site_id = 'site1', control = control)
# config <- getCloudConfig(site_id = 'site1')
# pdaSync(config)
```

Step 2: Derivative calculation

```r
# ############################  STEP 2: derivative  ###############################
## remote site3: please review and enter "1" to allow putting your derivatives to the server
i <- 3
# Sys.setenv(PDA_USER = 'site3', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
control <- pda(ipdata = lung_split[[i]], site_id = sites[i])
# now the group would see site3_derive.json	at https://pda.one

## remote site2: please review and enter "1" to allow putting your derivatives to the server
i <- 2
# Sys.setenv(PDA_USER = 'site2', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
control <- pda(ipdata = lung_split[[i]], site_id = sites[i])
# now the group would see site2_derive.json	at https://pda.one

## lead site1: please review and enter "1" to allow putting your derivatives to the server
i <- 1
# Sys.setenv(PDA_USER = 'site1', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
control <- pda(ipdata = lung_split[[i]], site_id = sites[i])
# now the group would see site1_derive.json	at https://pda.one
```

Step 3: Estimation

```r
# ############################  STEP 3: estimate  ###############################
## lead site1:
## please review and enter "1" to allow putting the surrogate estimate to the server
i <- 1
# Sys.setenv(PDA_USER = 'site1', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
control <- pda(ipdata = lung_split[[i]], site_id = sites[i])
# now the group would see site1_estimate.json	at https://pda.one
## compare the surrogate estimate with the pooled estimate
config <- getCloudConfig(site_id = 'site1')
fit.odal <- pdaGet(name = 'site1_estimate', config = config)
cbind(b.pool=fit.pool$coef,
      b.odal=fit.odal$btilde,
      sd.pool=summary(fit.pool)$coef[,2],
      sd.odal=sqrt(diag(solve(fit.odal$Htilde)/nrow(lung2))))
## the PDA ODAL is now completed!
## All the sites can still run their surrogate estimates and broadcast them.

## remote site2: (optional)
i <- 2
# Sys.setenv(PDA_USER = 'site2', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
control <- pda(ipdata = lung_split[[i]], site_id = sites[i])

## remote site3: (optional)
i <- 3
# Sys.setenv(PDA_USER = 'site3', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
control <- pda(ipdata = lung_split[[i]], site_id = sites[i])


## If all the sites broadcast their surrogate estimates,
## a final synthesize step can further improve the estimate.
## lead site1: uncomment to synchoronize the control before STEP 4
# Sys.setenv(PDA_USER = 'susernameite1', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
pda(site_id = 'site1', control = control)
config <- getCloudConfig(site_id = 'site1')
pdaSync(config)
```


Step 4: Synthesis (optional)

```r
# ########################  STEP 4: synthesize (optional)  ########################
## lead site1:
i <- 1
# Sys.setenv(PDA_USER = 'username', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
control <- pda(ipdata = lung_split[[i]], site_id = sites[i])
```



## References
1. Duan, R., Boland, M.R., Moore, J.H. and Chen, Y., (2019). [ODAL: a one-shot distributed algorithm to perform logistic regressions on electronic health records data from multiple clinical sites.](https://psb.stanford.edu/psb-online/proceedings/psb19/duan.pdf) *Pacific Symposium on Biocomputing* 2019 (pp. 30-41).
2. Duan, R., Boland, M., Liu, Z., Liu, Y., Chang, H., Xu., H, Chu, H., Schmid, C., Forrest, C., Holmes, J., Schuemie, M.J., Berlin, J.A., Moore, J.H. and Chen,Y., (2019). [Learning from electronic health records across multiple sites: a computationally and statistically efficient distributed algorithm.](https://pubmed.ncbi.nlm.nih.gov/31816040/) *Journal of the American Medical Informatics Association* 27(3), pp.376-385.
3. Duan, R., Luo, C., Schuemie, M.J., Tong, J., Liang, J., Boland, M.R., Bian, J., Xu, H., Berlin, J.A., Moore, J.H., Mahoney, K.B. and Chen, Y., (2020). [Learning from local to global - an efficient distributed algorithm for modeling time to event data.](https://doi.org/10.1093/jamia/ocaa044) *Journal of the American Medical Informatics Association* 27(7), pp.1028–1036
4. Tong, J., Duan, R., Li, R., Scheuemie, M.J., Moore, J.H. and Chen, Y., 2020, January. [Robust-ODAL: Learning from heterogeneous health systems without sharing patient-level data.](https://www.worldscientific.com/doi/abs/10.1142/9789811215636_0061) *In Pacific Symposium on Biocomputing. Pacific Symposium on Biocomputing* (Vol. 25, p. 695). NIH Public Access.
5. Duan, R., Chen, Z., Tong, J., Luo, C., Lyu, T., Tao, C., Maraganore, D., Bian, J. and Chen, Y., 2020. [Leverage Real-world Longitudinal Data in Large Clinical Research Networks for Alzheimer's Disease and Related Dementia (ADRD)](https://www.medrxiv.org/content/10.1101/2020.08.03.20167619v1). *medRxiv*.
6. Tong, J., Chen, Z., Duan, R., Lo-Ciganic, W.H., Lyu, T., Tao, C., Merkel, P.A., Kranzler, H.R., Bian, J. and Chen, Y., 2020. [Identifying Clinical Risk Factors for Opioid Use Disorder using a Distributed Algorithm to Combine Real-World Data from a Large Clinical Data Research Network](https://www.medrxiv.org/content/10.1101/2020.08.16.20167577v1). *medRxiv*.


## FAQ
#### *What do I need to do to run the package?*
#### *How to solve the errors about c++ in Rstudio?*
