# Copyright 2021 Penn Computing Inference Learning (PennCIL) lab
#       https://penncil.med.upenn.edu/team/
# This file is part of pda
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# https://style.tidyverse.org/functions.html#naming
# https://ohdsi.github.io/Hades/codeStyle.html#OHDSI_code_style_for_R


#' @useDynLib pda
#' @title Function to upload object to cloud as json
#' 
#' @usage pdaPut(obj,name,config,upload_without_confirm)
#' @param obj R object to encode as json and uploaded to cloud
#' @param name of file
#' @param config a list of variables for cloud configuration
#' @param upload_without_confirm logical. TRUE if want silent upload, no interactive confirm 
#' @importFrom utils menu
#' @return NONE
#' @seealso \code{pda}
#' @export
pdaPut <- function(obj,name,config,upload_without_confirm){
  obj_Json <- jsonlite::toJSON(obj, digits = 8)
  file_name <- paste0(name, '.json')
  
  if(!is.null(config$uri)){
    message(paste("Put",file_name,"on public cloud:"))
  }else{
    message(paste("Put",file_name,"on local directory", config$dir, ':'))
  }
  message(obj_Json)
  
  # if(interactive()) {
  if(upload_without_confirm==F) {
    authorize = menu(c("Yes", "No"), title="Allow transfer?")
  } else {
    authorize = "1"
  }
  if (authorize != 1) {
    warning("file not transferred. You can specify file transfer setting in pda().")
    return(FALSE)
  }
  # the file to upload
  if (is.character(config$dir)) {
    file_path <- paste0(config$dir,'/', file_name)
  } else {
    file_path <- paste0(tempdir(),'/', file_name)
  }
  write(obj_Json, file_path)
  if (is.character(config$uri)) {
    # create the url target of the file
    url <- file.path(config$uri, file_name)
    # webdav PUT request to send a file to cloud
    r<-httr::PUT(url, body = httr::upload_file(file_path), httr::authenticate(config$site_id, config$secret, 'digest'))
    message(paste("putting:",url))
  }
}

#' @useDynLib pda
#' @title Function to list available objects
#' @usage pdaList(config)
#' @param config a list of variables for cloud configuration
#' @importFrom rvest html_nodes
#' @import httr 
#' @return A list of (json) files on the cloud
#' @seealso \code{pda}
#' @export
pdaList <- function(config){
  if (is.character(config$uri)) {
    res<-httr::GET(config$uri, httr::authenticate(config$site_id, config$secret, 'digest'))
    files<-rvest::html_nodes(httr::content(res), xpath = "//table//td/a") 
    files<-files[lapply(files,length)>0]
    files<-regmatches(files, gregexpr("(?<=\")(.*?)(json)(?=\")", files, perl = TRUE))
  } else if (is.character(config$dir)) {
    files<-list.files(config$dir,pattern = "\\.json$") 
  } else {
    files<-list.files(tempdir(),pattern = "\\.json$") 
  }
  files<-substr(files,1,nchar(files)-5)
  return(files)
}


#' @useDynLib pda
#' @title Function to download json and return as object
#' @usage pdaGet(name,config)
#' @param name of file
#' @param config cloud configuration
#' @return A list of data objects from the json file on the cloud
#' @seealso \code{pda}
#' @export
pdaGet <- function(name,config){
  file_name <- paste0(name, '.json')
  #print(paste("Get",file_name,"from public cloud:"))
  # the file to upload
  if (is.character(config$dir)) {
    file_path <- paste0(config$dir,'/', file_name)
  } else {
    file_path <- paste0(tempdir(),'/', file_name)
  }
  if (is.character(config$uri)) {
    url <- file.path(config$uri, file_name)
    #write the file from GET request to file_path
    res<-httr::GET(url, httr::write_disk(file_path, overwrite = TRUE), httr::authenticate(config$site_id, config$secret, 'digest'))
    #print(paste("getting:",url))
  } 
  obj<-jsonlite::fromJSON(file_path)
  return(obj)
}


#' @useDynLib pda
#' @title gather cloud settings into a list
#' @usage getCloudConfig(site_id,dir,uri,secret)
#' @param site_id site identifier
#' @param dir shared directory path if flat files
#' @param uri web uri if web service
#' @param secret web token if web service
#' @return A list of cloud parameters: site_id, secret and uri
#' @seealso \code{pda}
#' @export
getCloudConfig <- function(site_id,dir=NULL,uri=NULL,secret=NULL){
  config<-list()
  pda_user<-Sys.getenv('PDA_USER')
  pda_secret<-Sys.getenv('PDA_SECRET')
  pda_uri<-Sys.getenv('PDA_URI')
  pda_dir<-Sys.getenv('PDA_DIR')
  config$site_id=site_id
  if(!is.null(secret)) {
    config$secret = secret
  } else if (pda_secret!='') {
    config$secret = pda_secret
  }
  
  if(!is.null(uri)) {
    config$uri = uri
  } else if (pda_uri!='') {
    config$uri = pda_uri
  } else{
    message('no cloud uri found! ')
  }
  
  if(!is.null(dir)) {
    config$dir = dir
  } else if (pda_dir!='') {
    config$dir = pda_dir
  }else{
    message('no public or local directory supplied, use local temporary:', tempdir())
    config$dir = tempdir()
  }
  config;
}



#' @useDynLib pda
#' @title PDA: Privacy-preserving Distributed Algorithm
#' 
#' @description  Fit Privacy-preserving Distributed Algorithms for linear, logistic, 
#'                Poisson and Cox PH regression with possible heterogeneous data across sites.
#' @usage pda(ipdata,site_id,control,dir,uri,secret)
#' @param ipdata  Local IPD data in data frame, should include at least one column for the outcome and one column for the covariates 
#' @param site_id Character site name
#' @param control pda control data
#' @param dir directory for shared flat file cloud
#' @param uri Universal Resource Identifier for this run
#' @param secret password to authenticate as site_id on uri
#' @param upload_without_confirm logical. TRUE if want silent upload, no interactive confirm 
#' @param hosdata (for dGEM) hospital-level data, should include the same name as defined in the control file
#' @return control
#' @seealso \code{pdaPut}, \code{pdaList}, \code{pdaGet}, \code{getCloudConfig} and \code{pdaSync}.
#' @import stats survival rvest jsonlite data.table httr Rcpp metafor
#'          
#' @references
#' Michael I. Jordan, Jason D. Lee & Yun Yang (2019) Communication-Efficient Distributed Statistical Inference, \cr
#'  \emph{Journal of the American Statistical Association}, 114:526, 668-681 \cr 
#'  \doi{10.1080/01621459.2018.1429274}.\cr 
#' (DLM) Yixin Chen, et al. (2006) Regression cubes with lossless compression and aggregation. 
#'    IEEE Transactions on Knowledge and Data Engineering, 18(12), pp.1585-1599. \cr
#' (DLMM) Chongliang Luo, et al. (2020) Lossless Distributed Linear Mixed Model with Application to Integration of Heterogeneous Healthcare Data.  
#'    medRxiv, \doi{10.1101/2020.11.16.20230730}. \cr
#' (DPQL) Chongliang Luo, et al. (2021) dPQL: a lossless distributed algorithm for generalized linear mixed model with application to privacy-preserving hospital profiling. \cr
#'    medRxiv, \doi{10.1101/2021.05.03.21256561}. \cr
#' (ODAL) Rui Duan, et al. (2020) Learning from electronic health records across multiple sites: \cr 
#'  A communication-efficient and privacy-preserving distributed algorithm. \cr 
#'  \emph{Journal of the American Medical Informatics Association}, 27.3:376–385,
#'  \cr \doi{10.1093/jamia/ocz199}.\cr 
#' (ODAC) Rui Duan, et al. (2020) Learning from local to global: An efficient distributed algorithm for modeling time-to-event data. \cr
#'   \emph{Journal of the American Medical Informatics Association}, 27.7:1028–1036, \cr 
#'    \doi{10.1093/jamia/ocaa044}. \cr
#' (ODACH) Chongliang Luo, et al. (2021) ODACH: A One-shot Distributed Algorithm for Cox model with Heterogeneous Multi-center Data. \cr
#'       \emph{medRxiv}, \doi{10.1101/2021.04.18.21255694}. \cr 
#' (ODAH) Mackenzie J. Edmondson, et al. (2021) An Efficient and Accurate Distributed Learning Algorithm for Modeling Multi-Site Zero-Inflated Count Outcomes. 
#'    medRxiv, pp.2020-12. \cr
#'    \doi{10.1101/2020.12.17.20248194}. \cr
#' (ADAP) Xiaokang Liu, et al. (2021) ADAP: multisite learning with high-dimensional heterogeneous data via A Distributed Algorithm for Penalized regression. \cr
#' (dGEM) Jiayi Tong, et al. (2022) dGEM: Decentralized Generalized Linear Mixed Effects Model \cr
#' @examples
#' require(survival)
#' require(data.table)
#' require(pda)
#' data(lung)
#' 
#' ## In the toy example below we aim to analyze the association of lung status with 
#' ## age and sex using logistic regression, data(lung) from 'survival', we randomly 
#' ## assign to 3 sites: 'site1', 'site2', 'site3'. we demonstrate using PDA ODAL can 
#' ## obtain a surrogate estimator that is close to the pooled estimate. We run the 
#' ## example in local directory. In actual collaboration, account/password for pda server 
#' ## will be assigned to the sites at the server https://pda.one.
#' ## Each site can access via web browser to check the communication of the summary stats.
#' 
#' ## for more examples, see demo(ODAC) and demo(ODAP)
#' 
#' # Create 3 sites, split the lung data amongst them
#' sites = c('site1', 'site2', 'site3')
#' set.seed(42)
#' lung2 <- lung[,c('status', 'age', 'sex')]
#' lung2$sex <- lung2$sex - 1
#' lung2$status <- ifelse(lung2$status == 2, 1, 0)
#' lung_split <- split(lung2, sample(1:length(sites), nrow(lung), replace=TRUE))
#' ## fit logistic reg using pooled data
#' fit.pool <- glm(status ~ age + sex, family = 'binomial', data = lung2)
#' 
#' 
#' # ############################  STEP 1: initialize  ###############################
#' control <- list(project_name = 'Lung cancer study',
#'                 step = 'initialize',
#'                 sites = sites,
#'                 heterogeneity = FALSE,
#'                 model = 'ODAL',
#'                 family = 'binomial',
#'                 outcome = "status",
#'                 variables = c('age', 'sex'),
#'                 optim_maxit = 100,
#'                 lead_site = 'site1',
#'                 upload_date = as.character(Sys.time()) )
#' 
#' 
#' ## run the example in local directory:
#' ## specify your working directory, default is the tempdir
#' mydir <- tempdir()
#' ## assume lead site1: enter "1" to allow transferring the control file  
#' pda(site_id = 'site1', control = control, dir = mydir)
#' ## in actual collaboration, account/password for pda server will be assigned, thus:
#' \dontrun{pda(site_id = 'site1', control = control, uri = 'https://pda.one', secret='abc123')}
#' ## you can also set your environment variables, and no need to specify them in pda:
#' \dontrun{Sys.setenv(PDA_USER = 'site1', PDA_SECRET = 'abc123', PDA_URI = 'https://pda.one')}
#' \dontrun{pda(site_id = 'site1', control = control)}
#' 
#' ##' assume remote site3: enter "1" to allow tranferring your local estimate 
#' pda(site_id = 'site3', ipdata = lung_split[[3]], dir=mydir)
#' 
#' ##' assume remote site2: enter "1" to allow tranferring your local estimate  
#' pda(site_id = 'site2', ipdata = lung_split[[2]], dir=mydir)
#' 
#' ##' assume lead site1: enter "1" to allow tranferring your local estimate  
#' ##' control.json is also automatically updated
#' pda(site_id = 'site1', ipdata = lung_split[[1]], dir=mydir)
#' 
#' ##' if lead site1 initialized before other sites,
#' ##' lead site1: uncomment to sync the control before STEP 2
#' \dontrun{pda(site_id = 'site1', control = control)}
#' \dontrun{config <- getCloudConfig(site_id = 'site1')}
#' \dontrun{pdaSync(config)}
#' 
#' #' ############################'  STEP 2: derivative  ############################ 
#' ##' assume remote site3: enter "1" to allow tranferring your derivatives  
#' pda(site_id = 'site3', ipdata = lung_split[[3]], dir=mydir)
#' 
#' ##' assume remote site2: enter "1" to allow tranferring your derivatives  
#' pda(site_id = 'site2', ipdata = lung_split[[2]], dir=mydir)
#' 
#' ##' assume lead site1: enter "1" to allow tranferring your derivatives  
#' pda(site_id = 'site1', ipdata = lung_split[[1]], dir=mydir)
#' 
#' 
#' #' ############################'  STEP 3: estimate  ############################ 
#' ##' assume lead site1: enter "1" to allow tranferring the surrogate estimate  
#' pda(site_id = 'site1', ipdata = lung_split[[1]], dir=mydir)
#' 
#' ##' the PDA ODAL is now completed!
#' ##' All the sites can still run their own surrogate estimates and broadcast them.
#' 
#' ##' compare the surrogate estimate with the pooled estimate 
#' config <- getCloudConfig(site_id = 'site1', dir=mydir)
#' fit.odal <- pdaGet(name = 'site1_estimate', config = config)
#' cbind(b.pool=fit.pool$coef,
#'       b.odal=fit.odal$btilde,
#'       sd.pool=summary(fit.pool)$coef[,2],
#'       sd.odal=sqrt(diag(solve(fit.odal$Htilde)/nrow(lung2))))
#'       
#' ## see demo(ODAL) for more optional steps
#' 
#' @return control
#' @export
pda <- function(ipdata=NULL,site_id,control=NULL,dir=NULL,uri=NULL,secret=NULL,
                upload_without_confirm=F,  
                hosdata=NULL # for dGEM
                ){ 
  config <- getCloudConfig(site_id,dir,uri,secret)
  #add a control if one was provided
  if(!(is.null(control)) &&  config$site_id==control$lead_site) { # control$sites[1]
    pdaPut(obj=control,name='control',config=config,upload_without_confirm)
    return(control)    # ?
  }
  control = pdaGet('control',config)
  message('You are performing Privacy-preserving Distributed Algorithm (PDA, https://github.com/Penncil/pda): ')
  message('your site = ', config$site_id)  
  
  if(control$model=='ODAL'){
    ODAL.steps<-c('initialize','derive','estimate','synthesize')
    ODAL.family<-'binomial'
  }else if(control$model=='ADAP'){
    ADAP.steps<-c('initialize','derive','estimate')
    ADAP.family<-'lasso'
  }else if(control$model=='ODAP'){
    ODAP.steps<-c('initialize','derive','estimate','synthesize')
    # for ODAP need to specify family (poisson, ztpoisson, quasipoisson, ztquasipoisson) in control
    ODAP.family<-control$family  
  }else if(control$model=='ODAPB'){
    ODAPB.steps<-c('initialize','derive','estimate','synthesize')
    # for ODAPB need to specify family (poisson) in control
    ODAPB.family<-control$family  
  }else if(control$model=='ODAH'){
    ODAH.steps<-c('initialize','derive','estimate','synthesize')
    #   family = 'hurdle' in control
    ODAH.family<-control$family 
    if(ODAH.family!='hurdle'){
      warning("currently only support family='hurdle'" )
      ODAH.family<-'hurdle'
    }
  }else if(control$model=='ODAC'){ 
    if(control$heterogeneity==F){     # ODAC
      ODAC.steps<-c('initialize','deriveUWZ','derive','estimate','synthesize')
    }else{                            # ODACH with heterogeneous baseline hazards across sites 
      ODAC.steps<-c('initialize','derive', 'estimate','synthesize')
    }
    ODAC.family<-'cox'
  }else if(control$model=='ODACAT'){  # multi-category  
    ODACAT.steps <- c('initialize','derive','estimate','synthesize')
    ODACAT.family <- 'multicategory'
  }else if(control$model=='ODACATH'){ # added by Jessie & Ken on Feb 24, 2023
    ODACATH.steps <- c('initialize','derive','estimate','synthesize')
    ODACATH.family <- 'multicategory'
    if(control$heterogeneity==T){
      message("You specified control$heterogeneity = T, so you are using the hetero-version of ODACAT.")
    }
  }else if(control$model=='DLM'){
    DLM.steps<-c('initialize', 'estimate')
    DLM.family<-'gaussian'
    if(control$heterogeneity==T){
      if(is.null(control$heterogeneity_effect)){
        stop('You specified control$heterogeneity = T, please also specify control$heterogeneity_effect as either "fixed" or "random"! ')
      } else if(control$heterogeneity_effect!='fixed' & control$heterogeneity_effect!='random'){
        stop('You specified control$heterogeneity = T, please also specify control$heterogeneity_effect as either "fixed" or "random"! ')
      } 
      if(length(setdiff(control$variables_heterogeneity, c(control$variables, "Intercept")))!=0)
        stop('You specified control$heterogeneity = T, please also specify control$variables_heterogeneity as a SUBSET of "Intercept" and control$variables!')
      if(is.null(control$variables_heterogeneity)){
        message('You specified control$heterogeneity = T, but no control$variables_heterogeneity, use "Intercept" as default!')
        control$variables_heterogeneity <- 'Intercept'      
      }
    }
  }else if(control$model=='DPQL'){
    # control$maxround: prespecified number of rounds
    # "derive_1"   "estimate_1" "derive_2"   "estimate_2" "derive_3"   "estimate_3"
    DPQL.steps<-paste0(rep(c('derive', 'estimate'), control$maxround), '_', rep(1:control$maxround, each=2))   
    DPQL.family<-control$family  # can be any glm family...
    # if(control$heterogeneity==T){
    #   if(is.null(control$heterogeneity_effect)){  # glmm with fixed site-specific effects?
    #     stop('You specified control$heterogeneity = T, please also specify control$heterogeneity_effect as either "fixed" or "random"! ')
    #   } else if(control$heterogeneity_effect!='fixed' & control$heterogeneity_effect!='random'){
    #     stop('You specified control$heterogeneity = T, please also specify control$heterogeneity_effect as either "fixed" or "random"! ')
    #   } 
    if(length(setdiff(control$variables_heterogeneity, c(control$variables, "Intercept")))!=0)
      stop('Please specify control$variables_heterogeneity as a SUBSET of "Intercept" and control$variables!')
    # stop('You specified control$heterogeneity = T, please also specify control$variables_heterogeneity as a SUBSET of "Intercept" and control$variables!')
    if(is.null(control$variables_heterogeneity)){
      message('No control$variables_heterogeneity, use "Intercept" as default!')
      # message('You specified control$heterogeneity = T, but no control$variables_heterogeneity, use "Intercept" as default!')
      control$variables_heterogeneity <- 'Intercept'      
    }
    # }
  }else if(control$model == 'dGEM'){
    dGEM.steps<-c('initialize','derive','estimate','synthesize')
    dGEM.family<-'binomial'
    variables_site_level <- control$variables_site_level
  }else if(control$model == 'OLGLM'){
    OLGLM.steps<-c('initialize')
    OLGLM.family<-'binomial'
  }else if(control$model == 'OLGLMM'){
    OLGLMM.steps<-c('initialize')
    OLGLMM.family<-'binomial'
  }else if(control$model == 'ODACH_CC'){ 
    ODACH_CC.steps<-c('initialize','derive', 'estimate','synthesize')
    ODACH_CC.family<-'cox'
  }
  
  family = get(paste0(control$model,'.family'))
  ## prepare the ipdata: make factor of all categorical variables to make dummy design matrix,
  ## in case some X's are degenerate at some site, see model.matrix(contrasts=...)
  # xlev.contrasts <- control$xlev
  # for(ii in names(control$xlev)){
  #   ipdata[,ii] = factor(ipdata[,ii], levels = control$xlev[[ii]])
  #   xlev.contrasts[[ii]] = 'contr.treatment' # options("contrasts") default
  # }
  
  if (!is.null(ipdata)){
    n = nrow(ipdata)
    if(family=='hurdle'){           # count and zero parts for hurdle, Xcount first
      variables <- control$variables_hurdle_count
    }else{
      variables <- control$variables
    }
    formula <- as.formula(paste(control$outcome, paste(variables, collapse = "+"), sep = '~'))
    mf <- model.frame(formula, ipdata, xlev=control$variables_lev)
  } 
  
  # Data sanity check: columns with no variation, or missed categorical levels...
  # if detected, may need to revise the data at this site, or 
  #     exclude problematic variables from the protocol, or
  #     exclude the site
  svd_d = svd(model.matrix(formula, mf))$d
  if(sum(svd_d < 1e-10)>=1) warning(site_id, ': data degeneration detected!!! Please discuss with your collaborators!')
  
  ## this is used in model.matrix(contrasts=...)
  # if(options()$contrasts['unordered']=="contr.treatment") options(contrasts = c("contr.treatment", "contr.poly"))
  
  # create ipdata via model.matrix to make dummy variables for categorical covariates...
  # the resulted ipdata format will be used in later functions, i.e. ODAX etc
  if(control$model=='ODAC'){  
    ipdata = data.table::data.table(time=as.numeric(model.response(mf))[1:n], 
                                    status=as.numeric(model.response(mf))[-c(1:n)], 
                                    model.matrix(formula, mf)[,-1])
    control$risk_factor = colnames(ipdata)[-c(1:2)]
  }else if(control$model=='ODAL'){
    ipdata = data.table::data.table(status=as.numeric(model.response(mf)), 
                                    model.matrix(formula, mf))
    control$risk_factor = colnames(ipdata)[-1]
  }else if(control$model=='ADAP'){
    ipdata = data.table::data.table(status=as.numeric(model.response(mf)), 
                                    model.matrix(formula, mf))
    control$risk_factor = colnames(ipdata)[-1]
  }else if(control$model=='ODAH'){  # count and zero parts for hurdle
    X_count = data.table::data.table(model.matrix(formula, mf))
    # also make design X_zero
    formula <- as.formula(paste(control$outcome, paste(control$variables_hurdle_zero, collapse = "+"), sep = '~'))
    mf <- model.frame(formula, ipdata)
    X_zero = data.table::data.table(model.matrix(formula, mf))
    if(is.character(control$offset)){
      offset <- ipdata[,control$offset]
    }else{
      offset = rep(0, n) 
    }
    ipdata <- list(ipdata=ipdata, X_count=X_count, X_zero=X_zero, offset=offset)  
    control$risk_factor = c('Intercept', control$variables_hurdle_count, 'Intercept', control$variables_hurdle_zero)   
  }else if(control$model=='ODAP'){
    ipdata = data.table::data.table(outcome=as.numeric(model.response(mf)), 
                                    offset=ifelse(is.character(control$offset), ipdata[,control$offset], 0),
                                    model.matrix(formula, mf))
    
    control$risk_factor = colnames(ipdata)[-c(1:2)]
  }else if(control$model=='ODAPB'){
    ipdata = data.table::data.table(outcome=as.numeric(model.response(mf)),
                                    offset=ifelse(is.character(control$offset), ipdata[,control$offset], 0),
                                    model.matrix(formula, mf))
    control$risk_factor = colnames(ipdata)[-c(1:2)]
  }else if(control$model=='ODACAT'){
    ipdata = data.table::data.table(outcome=as.numeric(model.response(mf)),  ## multi-category y is 1:q
                                    model.matrix(formula, mf))[,-2] # remove the intercept column. ODACAT does not need that column. 
    control$risk_factor = colnames(ipdata)[-1]
  }else if(control$model=='ODACATH'){ # added by Jessie & Ken on Feb 24, 2023
    ipdata = data.table::data.table(outcome=as.numeric(model.response(mf)),  ## multi-category y is 1:q
                                    model.matrix(formula, mf))[,-2] # remove the intercept column. ODACATH does not need that column. 
    control$risk_factor = colnames(ipdata)[-1]
  }else if(control$model=='DLM'){
    ipdata = data.table::data.table(outcome=as.numeric(model.response(mf)), 
                                    model.matrix(formula, mf))
    control$risk_factor = colnames(ipdata)[-1] 
    control$risk_factor_heterogeneity = control$risk_factor[grepl(paste0(control$variables_heterogeneity, collapse='|'), control$risk_factor)]
  }else if(control$model=='DPQL'){
    ipdata = data.table::data.table(outcome=as.numeric(model.response(mf)), 
                                    model.matrix(formula, mf))
    control$risk_factor = colnames(ipdata)[-1]           # may induce more cols for dummy vars
    control$risk_factor_heterogeneity = control$risk_factor[grepl(paste0(control$variables_heterogeneity, collapse='|'), control$risk_factor)]
  }else if(control$model=='dGEM'){
    if (!is.null(ipdata)){
      ipdata = data.table::data.table(status=as.numeric(model.response(mf)), 
                                      model.matrix(formula, mf))
      control$risk_factor = colnames(ipdata)[-1]
    }
  }else if(control$model=='OLGLM'){
    if(control$step == "initialize"){
      ipdata = data.table::data.table(status=as.numeric(model.response(mf)), 
                                      model.matrix(formula, mf))
      control$risk_factor = colnames(ipdata)[-1]
    }else{
      control$risk_factor = colnames(ipdata)[-1]
    }
  } else if(control$model=='OLGLMM'){
    if (!is.null(ipdata)){
      if(control$step == "initialize"){
        ipdata = data.table::data.table(status=as.numeric(model.response(mf)), 
                                        model.matrix(formula, mf))
        control$risk_factor = colnames(ipdata)[-1]
      }else{
        control$risk_factor = colnames(ipdata)[-1]
      }
    }
  } else if(control$model=='ODACH_CC'){
    if (!is.null(ipdata)){
      ipdata = data.table::data.table(time=as.numeric(model.response(mf))[1:n], 
                                      status=as.numeric(model.response(mf))[-c(1:n)],
                                      subcohort = ipdata$subcohort,
                                      # sampling_weight = ipdata$sampling_weight,
                                      model.matrix(formula, mf)[,-1])
      control$risk_factor = colnames(ipdata)[-c(1:3)] 
    }
  }
  
  
  if(is.character(control$step)){
    step_function <- paste0(control$model,'.', gsub('[^[:alpha:]]', '',control$step)) # "derive_1" for dPQL
    if(control$model == 'dGEM'){
      if(!is.null(ipdata)){
        if(control$step == 'derive'){
          step_obj <- get(step_function)(ipdata, control, config, hosdata)
        }else if (control$step == 'synthesize'){
          if(config$site_id != control$lead_site){
            stop("Only lead site or coordinating center can perform the last step (i.e., synthesize)")
          }else{
            step_obj <- get(step_function)(control, config)
          }
        } else {
          step_obj <- get(step_function)(ipdata, control, config)
        }
      }else {
        if (control$step == 'synthesize'){
          step_obj <- get(step_function)(control, config)
        }
      }
    }else if(control$model == "OLGLMM"){
      if(is.null(ipdata) & (config$site_id == control$lead_site)){
        print("As the leading site, you are going to produce the final results")
      }else{
        step_obj <- get(step_function)(ipdata, control, config)
      }
    }
    else{
      step_obj <- get(step_function)(ipdata, control, config)
    }
    
    
    if(control$step=='estimate'){
      if(control$model=='DLM'){
        message("Congratulations, the PDA is completed! The result is guaranteed to be identical to the pooled analysis")
      }else{
        if(control$model=='dGEM'){
          message("Congratulations, this the final step: you are transfering the counterfactural event rate. The lead site or coordinating center will broadcast the final results")
        }else{
          message("Congratulations, the PDA is completed! You can continue broadcasting your surrogate estimate to further synthesize them.")
        }
      }
      
    } else if(control$step==paste0('estimate_', control$maxround)){  # dPQL
      message("Congratulations, the PDA is completed! The result is guaranteed to be identical to the pooled analysis")
    }
    
    if(!is.null(ipdata)){
      pdaPut(step_obj,paste0(config$site_id,'_',control$step),config,upload_without_confirm)
    }else{
      if(control$model == "dGEM"){
        if(control$step == "synthesize"){
          pdaPut(step_obj,paste0(config$site_id,'_',control$step),config,upload_without_confirm)
        }
      }
    }
    
    #sync needed?
    if(config$site_id==control$lead_site) {
      control<-pdaSync(config,upload_without_confirm)
    }
  }
  invisible(control)
}


#' @useDynLib pda
#' @title pda control synchronize 
#' 
#' @description  update pda control if ready (run by lead)
#' @usage pdaSync(config)
#' @param config cloud configuration
#' @return control
#' @seealso \code{pda}
#' @export  
pdaSync <- function(config,upload_without_confirm){  
  control = pdaGet('control',config)
  if(control$model=='ODAL'){
    ODAL.steps<-c('initialize','derive','estimate','synthesize')
    ODAL.family<-'binomial'
  }else if(control$model=='ADAP'){
    ADAP.steps<-c('initialize','derive','estimate')
    ADAP.family<-'lasso'
  }else if(control$model=='ODAP'){
    ODAP.steps<-c('initialize','derive','estimate','synthesize')
    # for ODAP need to specify family (poisson, ztpoisson, quasipoisson, ztquasipoisson) in control
    ODAP.family<-control$family  
  }else if(control$model=='ODAPB'){
    ODAPB.steps<-c('initialize','derive','estimate','synthesize')
    # for ODAP need to specify family (poisson, ztpoisson, quasipoisson, ztquasipoisson) in control
    ODAPB.family<-control$family  
  }else if(control$model=='ODAH'){
    ODAH.steps<-c('initialize','derive','estimate','synthesize')
    # for ODAH family = 'hurdle' in control
    ODAH.family<-'hurdle' # control$family  
  }else if(control$model=='ODAC'){
    if(control$heterogeneity==F){     # ODAC
      ODAC.steps<-c('initialize','deriveUWZ','derive','estimate','synthesize')
    }else{                            # ODACH with heterogeneous baseline hazards across sites 
      ODAC.steps<-c('initialize','derive', 'estimate','synthesize')
    }
    ODAC.family<-'cox'
  } else if(control$model=='ODACAT'){
    ODACAT.steps<-c('initialize','derive','estimate','synthesize')
    ODACAT.family<-'multicategory'
  } else if(control$model=='ODACATH'){ # added by Jessie & Ken on Feb 24, 2023
    ODACATH.steps<-c('initialize','derive','estimate','synthesize')
    ODACATH.family<-'multicategory'
  } else if(control$model=='DLM'){
    DLM.steps<-c('initialize','estimate')
    DLM.family<-'gaussian'
  } else if(control$model=='DPQL'){
    DPQL.steps<-paste0(rep(c('derive', 'estimate'), control$maxround), '_', rep(1:control$maxround, each=2)) 
    DPQL.family<-control$family
  } else if(control$model=='dGEM'){
    dGEM.steps<-c('initialize','derive','estimate','synthesize')
    dGEM.family<-'binomial'
  }else if(control$model == 'OLGLM'){
    OLGLM.steps<-c('initialize')
    OLGLM.family<-'binomial'
  }else if(control$model == 'OLGLMM'){
    OLGLMM.steps<-c('initialize')
    OLGLMM.family<-'binomial'
  }else if(control$model=='ODACH_CC'){  # ODACH with case-cohort design
    ODACH_CC.steps<-c('initialize','derive', 'estimate','synthesize') 
    ODACH_CC.family<-'cox'
  }
  
  files<-pdaList(config) 
  if(all(paste0(control$sites,"_",control$step) %in% files)){ # all init are ready
    if(control$step=="initialize"){
      if(control$lead_site %in% control$sites){
        init_i <- pdaGet(paste0(control$lead_site,'_initialize'),config)
      }
      if(control$model=='DLM'){
        # DLM does not need derivative, thus estimate after initialize...
      }else if(control$model=='ODAH'){
        bhat_zero <-init_i$bhat_zero_i
        vbhat_zero <- init_i$Vhat_zero_i
        bhat_count <-init_i$bhat_count_i
        vbhat_count <- init_i$Vhat_count_i
        for(site_i in control$sites){
          if(site_i!=control$lead_site){
            init_i <- pdaGet(paste0(site_i,'_initialize'),config)
            bhat_zero = rbind(bhat_zero, init_i$bhat_zero_i)
            vbhat_zero = rbind(vbhat_zero, init_i$Vhat_zero_i)
            bhat_count = rbind(bhat_count, init_i$bhat_count_i)
            vbhat_count = rbind(vbhat_count, init_i$Vhat_count_i)
          }
        }
        #estimate from meta-analysis
        bmeta_zero = apply(bhat_zero/vbhat_zero,2,function(x){sum(x, na.rm = TRUE)})/apply(1/vbhat_zero,2,function(x){sum(x, na.rm = TRUE)})
        vmeta_zero = 1/apply(1/vbhat_zero,2,function(x){sum(x, na.rm = TRUE)})
        bmeta_count = apply(bhat_count/vbhat_count,2,function(x){sum(x, na.rm = TRUE)})/apply(1/vbhat_count,2,function(x){sum(x, na.rm = TRUE)})
        vmeta_count = 1/apply(1/vbhat_count,2,function(x){sum(x, na.rm = TRUE)})
        res = list(bmeta_zero = bmeta_zero, vmeta_zero = vmeta_zero, 
                   bmeta_count = bmeta_count, vmeta_count = vmeta_count)
        message('meta analysis (inverse variance weighted average) result:')
        #print(res)
        control$beta_zero_init = bmeta_zero
        control$beta_count_init = bmeta_count
      }else if(control$model=='ADAP'){
        bhat <- init_i$bhat_i 
        # vhbat = rep(nrow(ipdata), ncol(ipdata)-1) stores each site's sample size
        vbhat <- init_i$Vhat_i
        for(site_i in control$sites){
          if(site_i!=control$lead_site){
            init_i <- pdaGet(paste0(site_i,'_initialize'),config)
            bhat = rbind(bhat, init_i$bhat_i) 
            vbhat = rbind(vbhat, init_i$Vhat_i)
          }
        }
        #estimate with weighted average  
        bmeta = apply(diag(vbhat[,1])%*%bhat,2,function(x){sum(x, na.rm = TRUE)})/sum(vbhat[,1], na.rm = TRUE)
        vmeta = NA
        res = list(bmeta = bmeta, vmeta = vmeta)
        message('sample size weighted average result:')
        #print(res)
        control$beta_init = bmeta
      }else if(control$model == "OLGLM"){
        
        K <- length(control$sites)
        if(control$heterogeneity == FALSE){
          read_AD <- pdaGet(paste0(control$sites[1],'_initialize'),config)
          print(read_AD$SXY)
          
          SXY <- read_AD$SXY
          Xtable <- as.data.frame(matrix(unlist(read_AD$Xtable), ncol = (length(control$variables) + 2)))
          colnames(Xtable) <- c("intercept", control$variables, "n")
          
          Xcat <- Xtable[,colnames(Xtable)!='n']
          Xcat <- as.matrix(Xcat)
          counts <- Xtable$n
          for(site_i in control$sites[-1]){
            KSiteAD <- pdaGet(paste0(site_i,'_initialize'),config)
            SXY <- SXY+KSiteAD$SXY
            
            Xtable <- as.data.frame(matrix(unlist(KSiteAD$Xtable), ncol = (length(control$variables) + 2)))
            colnames(Xtable) <- c("intercept", control$variables, "n")
            counts <- counts + Xtable$n
          }
        }else{
          read_AD <- pdaGet(paste0(control$sites[1],'_initialize'),config)
          
          SXY.intercept <- read_AD$SXY[1]
          SXY.cov <- read_AD$SXY[-1]
          Xtable <-as.data.frame(matrix(unlist(read_AD$Xtable), ncol = (length(control$variables) + 2)))
          colnames(Xtable) <- c("intercept", control$variables, "n")
          
          Xcat0 <- Xtable[,colnames(Xtable)!='n']
          Xcat <- Xcat0[,-1]
          counts <- Xtable$n
          
          for(site_i in control$sites[-1]){
            KSiteAD <- pdaGet(paste0(control$sites[1],'_initialize'),config)
            SXY.intercept <- c(SXY.intercept,KSiteAD$SXY[1])
            SXY.cov <- SXY.cov+KSiteAD$SXY[-1]
            
            
            Xtable <- as.data.frame(matrix(unlist(KSiteAD$Xtable), ncol = (length(control$variables) + 2)))
            colnames(Xtable) <- c("intercept", control$variables, "n")
            
            Xcat0 <- Xtable[,colnames(Xtable)!='n']
            Xcat_tmp <- Xcat0[,-1]
            Xcat <- rbind(Xcat,Xcat_tmp)
            counts <- c(counts, Xtable$n)
          }
          
          SXY <- c(SXY.intercept,SXY.cov)
          SiteID <- rep(1:K, each = nrow(Xcat0))
          new.siteID <- sapply(1:K,function(i) ifelse(SiteID==i,1,0))
          colnames(new.siteID) <- paste0("Site", 1:K)
          Xcat <- cbind(new.siteID,Xcat)
          Xcat <- as.matrix(Xcat)
        }
        
        logLik_AD <- function(beta){
          -(sum(SXY*beta)-sum(log(1+exp(Xcat%*%c(beta)))*counts))/sum(counts)
        }
        
        fit.AD <- optim(par = rep(0, ncol(Xcat)), logLik_AD, method = "BFGS")
        
        se <- sqrt(diag(solve(hessian(func = function(x) logLik_AD(x)*sum(counts), x = fit.AD$par))))
        res <- data.frame(est = fit.AD$par, se = se)
        rownames(res) <- colnames(Xcat)
        
        control$final_output = res
        
      }else if(control$model == "OLGLMM"){
        for(site_i in control$site[1]){
          AD = pdaGet(paste0(site_i,'_initialize'),config)
          Xtable <- AD$Xtable
          SY <- Xtable$SY
          count <- Xtable$n
          new_Xtable <- Xtable[,1:(1+length(control$variables))]
          
          new_X <- new_Xtable[rep(seq_len(nrow(new_Xtable)), times = count), ]
          generate_Y <- function(Y_count, overall_count){
            sub_Y <- rep(c(1,0), times = c(Y_count, overall_count- Y_count))
          }
          new_Y <- unlist(mapply(generate_Y,SY,overall_count = count))
          
          output_0 = cbind(new_Y, new_X)
          colnames(output_0) = c(control$outcome,"intercept",control$variables)
          output_0$site <- site_i
        }
        
        for(site_i in control$sites[-1]){
          AD = pdaGet(paste0(site_i,'_initialize'),config)
          Xtable <- AD$Xtable
          SY <- Xtable$SY
          count <- Xtable$n
          new_Xtable <- Xtable[,1:(1+length(control$variables))]
          
          new_X <- new_Xtable[rep(seq_len(nrow(new_Xtable)), times = count), ]
          generate_Y <- function(Y_count, overall_count){
            sub_Y <- rep(c(1,0), times = c(Y_count, overall_count- Y_count))
          }
          new_Y <- unlist(mapply(generate_Y,SY,overall_count = count))
          
          output = cbind(new_Y, new_X)
          colnames(output) = c(control$outcome,"intercept",control$variables)
          output$site <- site_i
          output_0 <- rbind(output_0, output)
        }
        
        
        fit.pool.recons <- MASS::glmmPQL(as.formula(paste(control$outcome, paste(control$variables, collapse = "+"), sep = '~')), 
                                   ~1|site, 
                                   data = output_0,
                                   family='binomial')
        control$est <- fit.pool.recons$coefficients
        control$varFix <- fit.pool.recons$varFix
        control$rand_se <- fit.pool.recons$sigma
        
      }else{
        if(control$lead_site %in% control$sites){
          bhat <-init_i$bhat_i 
          vbhat <- init_i$Vhat_i
          if(control$model == "ODACATH"){
            bhat_eta = init_i$bhat_eta_i
          }
          for(site_i in control$sites){
            if(site_i!=control$lead_site){
              init_i <- pdaGet(paste0(site_i,'_initialize'),config)
              bhat = rbind(bhat, init_i$bhat_i)
              vbhat = rbind(vbhat, init_i$Vhat_i)
              if (control$model == "ODACATH"){
                bhat_eta = rbind(bhat_eta, init_i$bhat_eta_i)
              }
            }
          }
        }else{
          init_i = pdaGet(paste0(control$sites[1],'_initialize'),config)
          bhat <-init_i$bhat_i 
          vbhat <- init_i$Vhat_i
          
          for(site_i in control$sites[-1]){
            init_i <- pdaGet(paste0(site_i,'_initialize'),config)
            bhat = rbind(bhat, init_i$bhat_i)
            vbhat = rbind(vbhat, init_i$Vhat_i)
          }
        }
        #estimate from meta-analysis
        bmeta = apply(bhat/vbhat,2,function(x){sum(x, na.rm = TRUE)})/apply(1/vbhat,2,function(x){sum(x, na.rm = TRUE)})
        vmeta = 1/apply(1/vbhat,2,function(x){sum(x, na.rm = TRUE)})
        res = list(bmeta = bmeta, vmeta = vmeta)
        message('meta analysis (inverse variance weighted average) result:')
        
        ## sanity check: use more robust weighted median as init? 
        # bmeta = apply(bhat, 2, function(a) spatstat::weighted.median(a, site_size) )
        
        #print(res)
        control$beta_init = bmeta
        if (control$model == "ODACATH"){
          control$bhat_eta = bhat_eta
        }
        
      }
      mes <- 'beta_init added, step=2 (derivatives)! \n'
    }
    
    
    if(control$step=='derive'){
      if(control$model == "dGEM"){
        # get b_meta as initial bbar
        ghat <- c()
        vghat <- c()
        hosdata <- c()
        for(site_i in control$sites){
          i = 1
          init_i <- pdaGet(paste0(site_i,'_derive'),config)
          ghat = rbind(ghat, init_i$gammahat_i)
          vghat = rbind(vghat, init_i$Vgammahat_i)
          hosdata = rbind(hosdata, init_i$hosdata)
          i = i + 1
        }
        
        # meta-regression
        colnames(hosdata) = control$variables_site_level
        formula <- as.formula(paste("", paste(control$variables_site_level, collapse = "+"), sep = '~'))
        gamma_meta_reg_new = rma.uni(ghat, vghat, mods = formula, data = hosdata)
        gamma_BLUP <- blup(gamma_meta_reg_new)$pred
        
        control$estimated_hospital_effect = gamma_BLUP
      }else if(control$model == "OLGLM"){
        message("You are done!")
      }
    }
    
    
  }
  
  if(control$step=='synthesize'){
    if(control$model == "dGEM"){
      control$final_event_rate = pdaGet(paste0(control$lead_site,'_synthesize'),config)$final_event_rate
    }
  }
  # for DPQL, attach intermediate estimate to control after each round
  if(control$model=='DPQL' & control$step==paste0('estimate_', control$round)){
    est <- pdaGet(paste0(config$site_id,'_estimate_',control$round),config)
    control$bhat <- est$bhat
    control$uhat <- est$uhat
    if(control$round<control$maxround) control$round <- control$round + 1
  }
  
  steps = get(paste0(control$model,'.steps'))
  current_index <-  which(steps==control$step)
  if(current_index < length(steps)) {
    next_index <- current_index + 1
    next_step <- steps[next_index]
    mes <- paste0('next step=',next_index,' (',next_step,')! \n')
    control$step = next_step
  } else {
    mes <- paste0('finished! \n')
    control$step = NULL
  }
  message(mes)
  pdaPut(control,'control',config,upload_without_confirm)
  control
  
}
