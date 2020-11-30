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
- For ODAC (One-shot distributed algorithm for Cox regression), make sure you have cpp compiler as ODAC requires [Rcpp](https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-FAQ.pdf).


## Instructions for Installing and Running pda Package

Below are the instructions for installing and then running the package.

### How to install the pda package
There are several ways in which one could install the `pda` package. 

1. In RStudio, create a new project: File -> New Project... -> New Directory -> New Project. 

2. Execute the following R code: 

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

### How to run pda examples

Below are two ways to run the pda examples. 

In the toy example below we aim to analyze the association of lung status with age and sex using logistic regression, data(lung) from 'survival', we randomly assign to 3 sites: 'site1', 'site2', 'site3'. We demonstrate using PDA ODAL can obtain a surrogate estimator that is close to the pooled estimate. We run the example in local directory. In actual collaboration, account/password for pda server will be assigned to the sites at the server https://pda.one. Each site can access via web browser to check the communication of the summary stats.

You can either 

#### *Run example with demo()*

```r
demo(ODAL)
``` 
or
####  *Run example with code*

Step 0: load related R packages and prepare sample data

```r
# load packages
require(survival)
require(data.table)
require(pda)

# sample data, lung, from "survival" package
data(lung)

# create 3 sites, split the lung data amongst them
sites = c('site1', 'site2', 'site3')
set.seed(42)
lung2 <- lung[,c('time', 'status', 'age', 'sex')]
lung2$sex <- lung2$sex-1
lung2$status <- ifelse(lung2$status == 2, 1, 0)
lung_split<-split(lung2, sample(1:length(sites), nrow(lung), replace=TRUE))

## fit logistic reg using pooled data
fit.pool <- glm(status ~ age + sex, family = 'binomial', data = lung2)

``` 
Step 1: Initialization

```r
# ############################  STEP 1: initialize  ###############################
## lead site1: please review and enter "1" to allow putting the control file to the server
control <- list(project_name = 'Lung cancer study',
                step = 'initialize',
                sites = sites,
                heterogeneity = FALSE,
                model = 'ODAL',
                family = 'binomial',
                outcome = "status",
                variables = c('age', 'sex'),
                optim_maxit = 100,
                lead_site = 'site1',
                upload_date = as.character(Sys.time()) )

## run the example in local directory:
## assume lead site1: enter "1" to allow transferring the control file
pda(site_id = 'site1', control = control, dir = getwd())
## in actual collaboration, account/password for pda server will be assigned, thus:
# pda(site_id = 'site1', control = control, uri = 'https://pda.one', secret='abc123')
## you can also set your environment variables, and no need to specify them in pda:
# Sys.setenv(PDA_USER = 'site1', PDA_SECRET = 'abc123', PDA_URI = 'https://pda.one')
# pda(site_id = 'site1', control = control)

##' assume remote site3: enter "1" to allow tranferring your local estimate
pda(site_id = 'site3', ipdata = lung_split[[3]], dir=getwd())

##' assume remote site2: enter "1" to allow tranferring your local estimate
pda(site_id = 'site2', ipdata = lung_split[[2]], dir=getwd())


##' assume lead site1: enter "1" to allow tranferring your local estimate
##' control.json is also automatically updated
pda(site_id = 'site1', ipdata = lung_split[[1]], dir=getwd())

##' if lead site1 initialized before other sites,
##' lead site1: uncomment to sync the control before STEP 2
#' pda(site_id = 'site1', control = control)
#' config <- getCloudConfig(site_id = 'site1')
#' pdaSync(config)
``` 
Step 2: Calculate derivatives at each site

```r
#' ############################'  STEP 2: derivative  ###############################
##' assume remote site3: enter "1" to allow tranferring your derivatives
pda(site_id = 'site3', ipdata = lung_split[[3]], dir=getwd())

##' assume remote site2: enter "1" to allow tranferring your derivatives
pda(site_id = 'site2', ipdata = lung_split[[2]], dir=getwd())

##' assume lead site1: enter "1" to allow tranferring your derivatives
pda(site_id = 'site1', ipdata = lung_split[[1]], dir=getwd())
``` 
Step 3: Surrogate estimate

```r
#' ############################'  STEP 3: estimate  ###############################
##' assume lead site1: enter "1" to allow tranferring the surrogate estimate
pda(site_id = 'site1', ipdata = lung_split[[1]], dir=getwd())

##' the PDA ODAL is now completed!
##' All the sites can still run their own surrogate estimates and broadcast them.

``` 
Compare with the pooled and meta estimators

```r
##' compare the surrogate estimate with the pooled estimate
config <- getCloudConfig(site_id = 'site1', dir=getwd())
fit.odal <- pdaGet(name = 'site1_estimate', config = config)
control <- pdaGet(name = 'control', config)
cbind(b.pool=fit.pool$coef,
	   b.meta=control$beta_init,
      b.odal=fit.odal$btilde )

```

For other examples, please see 

```r
demo(ODAC)
```
for Cox regression, and 
 
```r
demo(ODAP)
```
for hurdle regression.




