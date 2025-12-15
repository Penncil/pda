# Copyright 2021 Penn Computing Inference Learning (PennCIL) lab
#       https://penncil.med.upenn.edu/
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

# rootSolve, prodlim,

#' @useDynLib pda
#' @title Function to upload object to cloud as json
#' 
#' @usage pdaPut(obj,name,config,upload_without_confirm=F,silent_message=F,digits=4)
#' @param obj R object to encode as json and uploaded to cloud
#' @param name of file
#' @param config a list of variables for cloud configuration
#' @param upload_without_confirm logical. TRUE if want silent upload, no interactive confirm 
#' @param silent_message logical. TRUE if want to mute message
#' @param digits digits after decimal points in the output json files
#' @importFrom utils menu
#' @return NONE
#' @seealso \code{pda}
#' @export
pdaPut <- function(obj,name,config,upload_without_confirm=F,silent_message=F,digits=16){
  mymessage <- function(mes, silent=silent_message) if(silent==F)  message(mes)
  
  obj_Json <- jsonlite::toJSON(obj, digits = digits)  # RJSONIO::toJSON(tt) keep vec name?
  file_name <- paste0(name, '.json')  
  
  # if(!is.null(config$uri)){
  #   mymessage(paste("Put",file_name,"on public cloud:"))
  # }else{
  #   mymessage(paste("Put",file_name,"on local directory", config$dir, ':'))
  # }
  # mymessage(obj_Json)
  
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
    mymessage(paste("putting:",url))
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
#' @usage getCloudConfig(site_id,dir=NULL,uri=NULL,secret=NULL,silent_message=T)
#' @param site_id site identifier
#' @param dir shared directory path if flat files
#' @param uri web uri if web service
#' @param secret web token if web service
#' @param silent_message logical, if the message will be muted
#' @return A list of cloud parameters: site_id, secret and uri
#' @seealso \code{pda}
#' @export
getCloudConfig <- function(site_id,dir=NULL,uri=NULL,secret=NULL,silent_message=T){
  mymessage <- function(mes, silent=silent_message) if(silent==F)  message(mes)
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
    # mymessage('no cloud uri found! ')
  }
  
  if(!is.null(dir)) {
    config$dir = dir
  } else if (pda_dir!='') {
    config$dir = pda_dir
  }else{
    mymessage(paste0('no public or local directory supplied, use local temporary:', tempdir()))
    config$dir = tempdir()
  }
  
  config;
}


# install.packages("data.tree")

#' @useDynLib pda
#' @title use this function to guide end-users step-by-step to identify best pda models for their tasks, and set up control.  
#' @usage pdaCatalog(task=c('Regression', 'Survival', 'Trial_emulation', 'Causal_inference', 'Design_analysis', 'Clustering'), write_json_file_path=getwd(), optim_maxit,optim_method,init_method)
#' @param task user-specified task, c('Regression', 'Survival', 'Trial_emulation', 'Causal_inference', 'Design_analysis', 'Clustering'). If no specify, display all models
#' @param write_json_file_path directory path to write the control file to
#' @param optim_maxit option in the control file for the optimization in pda, default 100
#' @param optim_method option in the control file for the optimization in pda, default "BFGS"
#' @param init_method option in the control file for calculating the initial estimate in pda, default "meta"
#' @return pda control 
#' @seealso \code{pda}
#' @export
pdaCatalog <- function(task=c('Regression', 
                              'Survival', 
                              'Trial_emulation', 
                              'Causal_inference', 
                              'Design_analysis', 
                              'Clustering' ), # 'Transfer_learning'
                       write_json_file_path=getwd(), 
                       optim_maxit,
                       optim_method,
                       init_method){
  # a = Hmisc::list.tree(object, fill = " | ", attr.print = F, size = F, maxlen = 1)     
  # library(data.tree)
  # library(jsonlite)
  # library(dplyr)
  S=readline(prompt="Here is the Catalog for all pda models. \nPlease choose the best pda model based on your Task and specifications. \nType  <Return>   to continue : ")  
  
  # Tree structure Catalog: Task - heterogeneity - heterogeneity_effect(0=fixed 1=random) or other specifications...
  catalog <- data.frame(
    pathString = c( 
      "Task",
      "Task/Regression, by outcome type",
      "Task/Regression, by outcome type/Continuous: DLM",
      "Task/Regression, by outcome type/Continuous: DLM/heterogeneity=No: DLM_0",
      "Task/Regression, by outcome type/Continuous: DLM/heterogeneity=Yes",
      "Task/Regression, by outcome type/Continuous: DLM/heterogeneity=Yes/heterogeneity_effect=fixed: DLM_1_0",
      "Task/Regression, by outcome type/Continuous: DLM/heterogeneity=Yes/heterogeneity_effect=random: DLM_1_1",
      "Task/Regression, by outcome type/Dichotomous",
      "Task/Regression, by outcome type/Dichotomous/# covariates > 30: ADAP",
      "Task/Regression, by outcome type/Dichotomous/counterfactual prediction: dGEM", # not involve dPQL
      # All covariates are categorical
      # "data fully stratified" -> "need one-shot lossless"? (need explanation... guideline book) 
      "Task/Regression, by outcome type/Dichotomous/data fully stratified", # Bingyu 
      "Task/Regression, by outcome type/Dichotomous/data fully stratified/heterogeneity=No: COLA_0_b",
      "Task/Regression, by outcome type/Dichotomous/data fully stratified/heterogeneity=Yes",
      "Task/Regression, by outcome type/Dichotomous/data fully stratified/heterogeneity=Yes/heterogeneity_effect=fixed: COLA_1_0_b",
      "Task/Regression, by outcome type/Dichotomous/data fully stratified/heterogeneity=Yes/heterogeneity_effect=random: COLA_1_1_b",
      # "Task/Regression, by outcome type/Dichotomous/heterogeneity=No", # set ODAL-R as built-in: if hessian has outlier then do median 
      "Task/Regression, by outcome type/Dichotomous/heterogeneity=No/need odds ratio: ODAL", 
      "Task/Regression, by outcome type/Dichotomous/heterogeneity=No/need risk ratio: ODAPB", 
      "Task/Regression, by outcome type/Dichotomous/heterogeneity=Yes: DPQL_b", 
      "Task/Regression, by outcome type/Count/excessive 0s: ODAH", 
      # "Task/Regression, by outcome type/count/non-zero counts: ODAPT",  # + truncated Pois? 
      # "Task/Regression, by outcome type/Count/excessive variation: DPLR", # Chongliang  
      "Task/Regression, by outcome type/Count/data fully stratified",     # the same as dichotomous
      "Task/Regression, by outcome type/Count/data fully stratified/heterogeneity=No: COLA_0_c",
      "Task/Regression, by outcome type/Count/data fully stratified/heterogeneity=Yes",
      "Task/Regression, by outcome type/Count/data fully stratified/heterogeneity=Yes/heterogeneity_effect=fixed: COLA_1_0_c",
      "Task/Regression, by outcome type/Count/data fully stratified/heterogeneity=Yes/heterogeneity_effect=random: COLA_1_1_c",
      "Task/Regression, by outcome type/Count/heterogeneity=No: ODAP",
      "Task/Regression, by outcome type/Count/heterogeneity=Yes: DPQL_c",
      "Task/Regression, by outcome type/Multi-category", 
      "Task/Regression, by outcome type/Multi-category/heterogeneity=No: ODACAT_0",
      "Task/Regression, by outcome type/Multi-category/heterogeneity=Yes: ODACAT_1", # ODACATH
      "Task/Survival",
      "Task/Survival/heterogeneity=No: ODAC_0",  
      "Task/Survival/heterogeneity=Yes: ODAC_1",  # ODACH
      "Task/Survival/time-varying effects: ODACT", # ODACT-H? 
      "Task/Survival/competing risk",         # Dazheng  
      "Task/Survival/competing risk/heterogeneity=No: ODACoR_0", #  
      "Task/Survival/competing risk/heterogeneity=Yes: ODACoR_1", #  
      "Task/Trial_emulation",      # will add LATTE T2E?
      "Task/Trial_emulation/propensity score stratificaion: LATTE_1",       # Lu   
      "Task/Trial_emulation/propensity score inverse weighting: LATTE_2",
      "Task/Trial_emulation/propensity score overlap weighting: LATTE_3",
      "Task/Causal_inference/Disc2o", # Jie TBA 
      "Task/Design_analysis/DRAFT", # Design-informed Regression Algorithm for Federated-learning Toolbox # ODACH-CC
      "Task/Clustering/ODEM",  # One-shot EM (Yudong)  
      "Task/Clustering/ODMM",  # 
      "Task/Clustering/DMLCA"  # Distributed EM (Xiaokang)
      # "Task/Transfer_learning"  # Jie TBA
    ) # , value = LETTERS[1:16] # Optional: Add values to nodes
  ) 
  catalog_tree <- data.tree::as.Node(catalog)
  print(catalog_tree) # how to change "levelName" to 'pda model'?
  
  # read in the model name, as identified by the pda Catalog
  model = readline(prompt='\nPlease provide the pda model as identified by the Catalog, \nif specfication (e.g. _1_1 after DLM) is not provided, default will be used: ')
  
  ## Start generate pda control
  control = list()
  
  # read in model as selected from the Catalog
  model = unlist(strsplit(model, '_'))
  tt=as.character(model[2])
  control$model = model[1]
  control$heterogeneity = ifelse(model[2]==0, F, T)
  control$heterogeneity_effect = ifelse(model[3]==0, 'fixed', 'random')
  if(model[1]=='LATTE') control$propensity_score = dplyr::case_when(model[2]==1 ~ 'stratificaion',
                                                             model[2]==2 ~ 'inverse weighting',
                                                             model[2]==3 ~ 'overlap weighting')
  
  # family from model
  fml = dplyr::case_when(model[1]%in%c('DLM') ~ 'gaussian',
                  model[1]%in%c('ODAL','dGEM','ODAPB') ~ 'binomial',
                  model[1]%in%c('ODAP','DPLR') ~ 'poisson',
                  model[1]%in%c('ODACAT') ~ 'multicategory',
                  model[1]%in%c('ODAC','ODACT','ODACoR') ~ 'cox',
                  model[1]%in%c('ADAP') ~ 'lasso',
                  model[1]%in%c('ODAH') ~ 'hurdle')
  # for COLA: control$mixed_effects?
  if(model[1] %in%c('DPQL', 'COLA') ) fml = dplyr::case_when(model[length(model)]=='b' ~ 'binomial',
                                                      model[length(model)]=='c' ~ 'poisson' )
  control$family = fml
  
  # read in project_name
  control$project_name = readline(prompt="\nPlease name your project: ")
  
  # read in site names
  sites = readline(prompt='\nPlease list the participating sites separated by blanks, \nwith the first site as lead/coordinating site, \ntype <Return> to skip if you want to specify it later: ')
  control$sites = unlist(strsplit(sites, ' ', fixed = T))
  control$lead_site = control$sites[1]
  
  # read in variable names: 
  # outcome variable
  outcome = readline(prompt='\nPlease provide your outcome variable name, \nif your Task is Survival, provide both variables for time-to-event and censor status separated by blank, \ntype <Return> to skip if Task is Clustering: ') 
  if(grepl('ODAC', model[1])){
    outcome = unlist(strsplit(outcome, ' ', fixed = T))
    control$outcome = paste0('Surv(', outcome[1], ', ', outcome[2], ')' )
  } 
  
  # covariate variables
  variables = readline(prompt='\nPlease provide your covariate variable names separated by blanks, \nthe variables need to be harmonized properly between sites (see Harminization guidelines xxx), \ntype <Return> to skip if you want to provide it later: ')
  control$variables = unlist(strsplit(variables, ' ', fixed = T))
  
  # levels of all categorical X's, with the first being the reference
  variables_lev = readline(prompt="\nPlease provide the levels of all categorical covariates, with the first being the reference, \n example input: age=c('young', 'old'), sex=c('M','F') \ntype <Return> to skip if you want to provide it later: ")
  control$variables_lev = eval(parse(text=paste0('list(', variables_lev, ')')))
  
  # other defaults
  control$step = 'initialize'
  # control$family = ''  # from model?...
  # ask for input if user is statistician?
  control$optim_maxit = 100
  control$optim_method = "BFGS"
  control$init_method = 'meta' 
  control$upload_date = Sys.time()
  
  write_json_file = menu(c("Yes", "No"), title="\nDo you want to write control into a json file at your specified folder? ")  
  if(write_json_file==1) write(toJSON(control, digits = 4), file=paste0(write_json_file_path, '/control.json'))
  
  return(control) 
}


#' @useDynLib pda
#' @title PDA: Privacy-preserving Distributed Algorithm
#' 
#' @description  Fit Privacy-preserving Distributed Algorithms for linear, logistic, 
#'                Poisson and Cox PH regression with possible heterogeneous data across sites.
#' @usage pda(ipdata=NULL,site_id,control=NULL,dir=NULL,uri=NULL,secret=NULL,upload_without_confirm=F, silent_message=F, digits=4,hosdata=NULL)
#' @param ipdata  Local IPD data in data frame, should include at least one column for the outcome and one column for the covariates 
#' @param site_id Character site name
#' @param control pda control data
#' @param dir directory for shared flat file cloud
#' @param uri Universal Resource Identifier for this run
#' @param secret password to authenticate as site_id on uri
#' @param upload_without_confirm logical. TRUE if want silent upload, no interactive confirm 
#' @param silent_message logical. TRUE if want to mute message
#' @param digits digits after decimal points in the output json files
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
#' (COLA) Wu, Q., Reps, J.M., Li, L. et al. COLA-GLM: collaborative one-shot and lossless algorithms of generalized linear models for decentralized observational healthcare data. npj Digit. Med. 8, 442 (2025). https://doi.org/10.1038/s41746-025-01781-1. \cr
#' (ODACT) Liang CJ, Luo C, Kranzler HR, Bian J, Chen Y. Communication-efficient federated learning of temporal effects on opioid use disorder with data from distributed research networks. J Am Med Inform Assoc. 2025 Apr 1;32(4):656-664. doi: 10.1093/jamia/ocae313. PMID: 39864407; PMCID: PMC12005629. \cr
#' (DisC2o) Tong J, et al. 2025. DisC2o-HD: Distributed causal inference with covariates shift for analyzing real-world high-dimensional data. Journal of Machine Learning Research. 2025;26(3):1-50. \cr
#' @return control
#' @export
pda <- function(ipdata=NULL,site_id,control=NULL,dir=NULL,uri=NULL,secret=NULL,
                upload_without_confirm=F, silent_message=F, digits=16,
                hosdata=NULL # for dGEM
){ 
  config <- getCloudConfig(site_id,dir,uri,secret,silent_message)
  mymessage <- function(mes, silent=silent_message) if(silent==F)  message(mes)
  files <- pdaList(config)
  # mymessage('You are performing Privacy-preserving Distributed Algorithm (PDA, https://github.com/Penncil/pda): ')
  mymessage(paste0('your site = ', config$site_id)) 
  
  data_now<-ipdata
  # read in control, or lead site add a control file to the cloud if there is none
  if('control' %in% files) {
    control = pdaGet('control',config) 
  } else { 
    if(!(is.null(control)) &&  config$site_id==control$lead_site) {    
      pdaPut(obj=control,name='control',config=config,upload_without_confirm,silent_message,digits)
      return(control)
    } else {
      stop('A control file is needed from the lead site!') 
    }
  } 
  
  ## specify pda steps and family based on model
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
      mymessage("You specified control$heterogeneity = T, so you are using the hetero-version of ODACAT.")
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
        mymessage('You specified control$heterogeneity = T, but no control$variables_heterogeneity, use "Intercept" as default!')
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
      mymessage('No control$variables_heterogeneity, use "Intercept" as default!')
      # message('You specified control$heterogeneity = T, but no control$variables_heterogeneity, use "Intercept" as default!')
      control$variables_heterogeneity <- 'Intercept'      
    } 
  }else if(control$model == 'dGEM'){
    dGEM.steps<-c('initialize','derive','estimate','synthesize')
    dGEM.family<-'binomial'
    variables_site_level <- control$variables_site_level
  }else if(control$model == 'COLA'){
    COLA.steps <- c('initialize', 'estimate')
    COLA.family <- control$family
    if (control$mixed_effects==TRUE){
      mymessage('You are setting the existence of mixed effects (mixed_effects = TRUE) and assuming to implement COLA-GLMM.')
    } 
  }else if(control$model == 'ODACH_CC'){ 
    ODACH_CC.steps <- c('initialize','derive', 'estimate','synthesize')
    ODACH_CC.family <- 'cox'
  }else if(control$model=='DisC2o'){
    DisC2o.steps <- c('PSinitialize','PSderive','PSestimate',
                      'OMinitialize','OMderive','OMestimate',
                      'AIPWestimate','synthesize')
    DisC2o.family <- control$family
  } else if(control$model=='ODACT'){ # ODACH with time-varying effects
    ODACT.steps <- c('initialize','derive', 'estimate','synthesize')
    ODACT.family <- 'cox'
  } else if(control$model=='LATTE'){  
    LATTE.steps<-c('initialize','estimate')
    LATTE.family<-'binomial'
  }
  
  family = get(paste0(control$model,'.family'))
  
  ## prepare the ipdata: make factor of all categorical variables to make dummy design matrix,
  # in case some X's are degenerate at some site, see model.matrix(contrasts=...)
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
    if(control$model=="DisC2o"){
      treatment_name <- control$treatment
      treatment_col <- ipdata[[treatment_name]]
    }
    formula <- as.formula(paste(control$outcome, paste(variables, collapse = "+"), sep = '~'))
    mf <- model.frame(formula, ipdata, xlev=control$variables_lev)
  } 
  
  # Data sanity check: columns with no variation, or missed categorical levels...
  # if detected, may need to revise the data at this site, or 
  #     exclude problematic variables from the protocol, or
  #     exclude the site
  svd_d = svd(model.matrix(formula, mf))$d
  if(sum(svd_d < 1e-10)>=1) warning(site_id, ': data degeneration detected!!! Proceed only if this is expected!')
  
  ## this is used in model.matrix(contrasts=...)
  # if(options()$contrasts['unordered']=="contr.treatment") options(contrasts = c("contr.treatment", "contr.poly"))
  
  # create ipdata via model.matrix to make dummy variables for categorical covariates...
  # the resulted ipdata format will be used in later functions, i.e. ODAX etc
  if(control$model=='ODAC'){  
    ipdata = data.table::data.table(time=as.numeric(model.response(mf))[1:n], 
                                    status=as.numeric(model.response(mf))[-c(1:n)], 
                                    model.matrix(formula, mf)[,-1])
    ipdata = data.table(data.frame(ipdata)) 
    control$risk_factor = colnames(ipdata)[-c(1:2)]
  }else if(control$model=='ODAL'){
    ipdata = data.table::data.table(status=as.numeric(model.response(mf)), 
                                    model.matrix(formula, mf))
    ipdata = data.table(data.frame(ipdata)) 
    control$risk_factor = colnames(ipdata)[-1]
  }else if(control$model=='ADAP'){
    ipdata = data.table::data.table(status=as.numeric(model.response(mf)), 
                                    model.matrix(formula, mf))
    ipdata = data.table(data.frame(ipdata)) 
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
    ipdata = data.table(data.frame(ipdata)) 
    control$risk_factor = colnames(ipdata)[-c(1:2)]
  }else if(control$model=='ODAPB'){
    ipdata = data.table::data.table(outcome=as.numeric(model.response(mf)),
                                    offset=ifelse(is.character(control$offset), ipdata[,control$offset], 0),
                                    model.matrix(formula, mf))
    ipdata = data.table(data.frame(ipdata)) 
    control$risk_factor = colnames(ipdata)[-c(1:2)]
  }else if(control$model=='ODACAT'){
    ipdata = data.table::data.table(outcome=as.numeric(model.response(mf)),  ## multi-category y is 1:q
                                    model.matrix(formula, mf))[,-2] # remove the intercept column. ODACAT does not need that column. 
    ipdata = data.table(data.frame(ipdata)) 
    control$risk_factor = colnames(ipdata)[-1]
  }else if(control$model=='ODACATH'){ # added by Jessie & Ken on Feb 24, 2023
    ipdata = data.table::data.table(outcome=as.numeric(model.response(mf)),  ## multi-category y is 1:q
                                    model.matrix(formula, mf))[,-2] # remove the intercept column. ODACATH does not need that column. 
    ipdata = data.table(data.frame(ipdata)) 
    control$risk_factor = colnames(ipdata)[-1]
  }else if(control$model=='DLM'){
    ipdata = data.table::data.table(outcome=as.numeric(model.response(mf)), 
                                    model.matrix(formula, mf))
    ipdata = data.table(data.frame(ipdata)) 
    control$risk_factor = colnames(ipdata)[-1] 
    control$risk_factor_heterogeneity = control$risk_factor[grepl(paste0(control$variables_heterogeneity, collapse='|'), control$risk_factor)]
  }else if(control$model=='DPQL'){
    ipdata = data.table::data.table(outcome=as.numeric(model.response(mf)), 
                                    model.matrix(formula, mf))
    ipdata = data.table(data.frame(ipdata)) 
    control$risk_factor = colnames(ipdata)[-1]           # may induce more cols for dummy vars
    control$risk_factor_heterogeneity = control$risk_factor[grepl(paste0(control$variables_heterogeneity, collapse='|'), control$risk_factor)]
  }else if(control$model=='dGEM'){
    if (!is.null(ipdata)){
      ipdata = data.table::data.table(status=as.numeric(model.response(mf)), 
                                      model.matrix(formula, mf))
      ipdata = data.table(data.frame(ipdata)) 
      control$risk_factor = colnames(ipdata)[-1]
    } 
  }else if (control$model == 'COLA'){
    if (isTRUE(control$mixed_effects)) {
      # --- COLA-GLMM ---
      keep <- c(control$variables, control$outcome)
      # missing <- setdiff(keep, names(ipdata))
      # if (length(missing)) stop("Missing columns for COLA-GLMM: ", paste(missing, collapse = ", "))
      
      # ensure data.table and select using ..keep
      ipdata <- data.table::as.data.table(ipdata)[, ..keep]
      
      # standardize outcome name for COLA.initialize() expectations
      setnames(ipdata, control$outcome, control$outcome, skip_absent = TRUE)
      
      # (optional) coerce outcome to numeric (either binary or counts)
      if (!is.numeric(ipdata[[control$outcome]])) {
        ipdata[,(control$outcome) := as.numeric(get(control$outcome))]
      }
      
      # Do NOT set control$risk_factor here; COLA.initialize builds x_names/X0 itself.
    } else {
      # --- COLA-GLM / COLA-GLM-H path for design-matrix ---
      if (control$step == "initialize") {
        ipdata <- data.table::data.table(
          outcome = as.numeric(model.response(mf)),
          model.matrix(formula, mf)
        )
        control$risk_factor <- colnames(ipdata)[-1]
      } else {
        control$risk_factor <- colnames(ipdata)[-1]
      }
    }
  }else if(control$model=='ODACH_CC'){
    if (!is.null(ipdata)){
      if(ncol(model.response(mf))==3){
        ipdata = data.table::data.table(time_in=as.numeric(model.response(mf)[,1]), 
                                        time=as.numeric(model.response(mf)[,2]), 
                                        status=as.numeric(model.response(mf)[,3]),
                                        subcohort = ipdata$subcohort,
                                        # sampling_weight = ipdata$sampling_weight,
                                        model.matrix(formula, mf)[,-1])
      } else if(ncol(model.response(mf))==2){
        ipdata = data.table::data.table(time=as.numeric(model.response(mf))[1:n], 
                                        status=as.numeric(model.response(mf))[-c(1:n)],
                                        subcohort = ipdata$subcohort,
                                        # sampling_weight = ipdata$sampling_weight,
                                        model.matrix(formula, mf)[,-1])
      }
      # convert irregular risk factor names, e.g. `Group (A,B,C) B` to Group..A.B.C..B
      # this should (and will) apply to all other models...
      ipdata = data.table(data.frame(ipdata)) 
      control$risk_factor = colnames(ipdata)[-c(1:(ncol(model.response(mf))+1))] 
    }
  }else if(control$model=='DisC2o'){
    ipdata = data.table::data.table(treatment = treatment_col,
                                    status=as.numeric(model.response(mf)), 
                                    model.matrix(formula, mf))
    ipdata = data.table(data.frame(ipdata)) 
    control$risk_factor = colnames(ipdata)[-c(1,2)]
  }else if(control$model=='ODACT'){  
    ipdata = data.table::data.table(time=as.numeric(model.response(mf))[1:n], 
                                    status=as.numeric(model.response(mf))[-c(1:n)], 
                                    model.matrix(formula, mf)[,-1])
    ipdata = data.table(data.frame(ipdata)) 
    control$risk_factor = colnames(ipdata)[-c(1:2)]
  }
  
  ## synchronize control file (at lead site), if lead site sees all collab sites ready
  files<-pdaList(config) 
  if(config$site_id==control$lead_site & all(paste0(control$sites,"_",control$step) %in% files)) {
    control<-pdaSync(config,upload_without_confirm,silent_message,digits)
  }
  
  ## execute the current step function
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
    }
    else{
      step_obj <- get(step_function)(ipdata, control, config)
    }
    
    ## pda completion message
    if(control$step=='estimate'){
      if(control$model=='DLM'){
        mymessage("Congratulations, the PDA is completed! The result is guaranteed to be identical to the pooled analysis")
      }else if(control$model=='LATTE'){
        if(config$site_id==control$lead_site) {
          mymessage("Congratulations, the PDA is completed! The result is guaranteed to be identical to the pooled analysis")
        }
    } else{
        if(control$model=='dGEM'){
          mymessage("Congratulations, this the final step: you are transfering the counterfactural event rate. The lead site or coordinating center will broadcast the final results")
        }else{
          mymessage("Congratulations, the PDA is completed! You can continue broadcasting your surrogate estimate to further synthesize them.")
        }
      }
    } else if(control$step==paste0('estimate_', control$maxround)){  # dPQL
      mymessage("Congratulations, the PDA is completed! The result is guaranteed to be identical to the pooled analysis")
    }
    
    ## write output to .json file
    if(!is.null(ipdata)){
      pdaPut(step_obj,paste0(config$site_id,'_',control$step),config,upload_without_confirm,silent_message,digits)
    }else{
      if(control$model == "dGEM"){
        if(control$step == "synthesize"){
          pdaPut(step_obj,paste0(config$site_id,'_',control$step),config,upload_without_confirm,silent_message,digits)
        }
      }
    }  
    
    if((control$step=="PSinitialize" & site_id==control$lead_site) |
       (control$step=="OMinitialize" & site_id==control$lead_site )) {
      control<-pdaSync(config,upload_without_confirm,silent_message,digits)
    }
    
    if(control$step=="PSestimate" | control$step=="OMestimate" ){
      control<-pdaSync(config,upload_without_confirm,silent_message,digits)
      pda(site_id = site_id, ipdata = data_now, dir=dir)
    }
    
    
    ## synchronize control file (at lead site), if lead site sees all collab sites ready
    files<-pdaList(config) 
    if(config$site_id==control$lead_site & all(paste0(control$sites,"_",control$step) %in% files)) {
      if(control$model=='DisC2o'& site_id==control$lead_site){
        control<-pdaSync(config,upload_without_confirm,silent_message,digits)
        pda(site_id = site_id, ipdata = data_now, dir=dir)
      }else{
        control<-pdaSync(config,upload_without_confirm,silent_message,digits)
      }
    }     
    
  }
  invisible(control)
}


#' @useDynLib pda
#' @title pda control synchronize 
#' 
#' @description  update pda control if ready (run by lead)
#' @usage pdaSync(config,upload_without_confirm,silent_message, digits)
#' @param config cloud configuration
#' @param upload_without_confirm logical. TRUE if want silent upload, no interactive confirm 
#' @param silent_message logical. TRUE if want to mute message
#' @param digits digits after decimal points in the output json files
#' @return control
#' @seealso \code{pda}
#' @export  
pdaSync <- function(config,upload_without_confirm,silent_message=F, digits=16){  
  control = pdaGet('control',config)
  mymessage <- function(mes, silent=silent_message) if(silent==F)  message(mes)
  
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
  }else if(control$model == 'COLA'){
    COLA.steps <- c('initialize','estimate')
    COLA.family <- control$family
  }else if(control$model=='ODACH_CC'){  # ODACH with case-cohort design
    ODACH_CC.steps <- c('initialize','derive', 'estimate','synthesize') 
    ODACH_CC.family <- 'cox'
  }else if (control$model == "LATTE") {
    LATTE.steps <- c("initialize", "estimate")
    LATTE.family <- "binomial"
  }else if(control$model=='DisC2o'){
    DisC2o.steps <- c('PSinitialize','PSderive','PSestimate',
                      'OMinitialize','OMderive','OMestimate',
                      'AIPWestimate','synthesize')
    DisC2o.family <- control$family
  }else if(control$model=='ODACT'){
    ODACT.steps <- c('initialize','derive', 'estimate','synthesize') 
    ODACT.family <- 'cox'
  }
  
  files<-pdaList(config) 
  if(all(paste0(control$sites,"_",control$step) %in% files)){ # all init are ready
    if(control$step=="initialize"){
      if(control$lead_site %in% control$sites){
        init_i <- pdaGet(paste0(control$lead_site,'_initialize'),config)
      }
      if(control$model=='DLM'){
        # DLM does not need derivative, thus estimate after initialize...
      }else if(control$model=='LATTE'){
        # LATTE does not need derivative, thus estimate after initialize...
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
        mymessage('meta analysis (inverse variance weighted average) result:')
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
        mymessage('sample size weighted average result:')
        #print(res)
        control$beta_init = bmeta
        # }else if(control$model == "OLGLM"){ 
        #   K <- length(control$sites)
        #   if(control$heterogeneity == FALSE){
        #     read_AD <- pdaGet(paste0(control$sites[1],'_initialize'),config)
        #     # print(read_AD$SXY) 
        #     SXY <- read_AD$SXY
        #     Xtable <- as.data.frame(matrix(unlist(read_AD$Xtable), ncol = (length(control$variables) + 2)))
        #     colnames(Xtable) <- c("intercept", control$variables, "n") 
        #     Xcat <- Xtable[,colnames(Xtable)!='n']
        #     Xcat <- as.matrix(Xcat)
        #     counts <- Xtable$n
        #     for(site_i in control$sites[-1]){
        #       KSiteAD <- pdaGet(paste0(site_i,'_initialize'),config)
        #       SXY <- SXY+KSiteAD$SXY 
        #       Xtable <- as.data.frame(matrix(unlist(KSiteAD$Xtable), ncol = (length(control$variables) + 2)))
        #       colnames(Xtable) <- c("intercept", control$variables, "n")
        #       counts <- counts + Xtable$n
        #     }
        #   }else{
        #     read_AD <- pdaGet(paste0(control$sites[1],'_initialize'),config) 
        #     SXY.intercept <- read_AD$SXY[1]
        #     SXY.cov <- read_AD$SXY[-1]
        #     Xtable <-as.data.frame(matrix(unlist(read_AD$Xtable), ncol = (length(control$variables) + 2)))
        #     colnames(Xtable) <- c("intercept", control$variables, "n") 
        #     Xcat0 <- Xtable[,colnames(Xtable)!='n']
        #     Xcat <- Xcat0[,-1]
        #     counts <- Xtable$n
        #     
        #     for(site_i in control$sites[-1]){
        #       KSiteAD <- pdaGet(paste0(control$sites[1],'_initialize'),config)
        #       SXY.intercept <- c(SXY.intercept,KSiteAD$SXY[1])
        #       SXY.cov <- SXY.cov+KSiteAD$SXY[-1] 
        #       Xtable <- as.data.frame(matrix(unlist(KSiteAD$Xtable), ncol = (length(control$variables) + 2)))
        #       colnames(Xtable) <- c("intercept", control$variables, "n") 
        #       Xcat0 <- Xtable[,colnames(Xtable)!='n']
        #       Xcat_tmp <- Xcat0[,-1]
        #       Xcat <- rbind(Xcat,Xcat_tmp)
        #       counts <- c(counts, Xtable$n)
        #     }
        #     
        #     SXY <- c(SXY.intercept,SXY.cov)
        #     SiteID <- rep(1:K, each = nrow(Xcat0))
        #     new.siteID <- sapply(1:K,function(i) ifelse(SiteID==i,1,0))
        #     colnames(new.siteID) <- paste0("Site", 1:K)
        #     Xcat <- cbind(new.siteID,Xcat)
        #     Xcat <- as.matrix(Xcat)
        #   } 
        #   logLik_AD <- function(beta){
        #     -(sum(SXY*beta)-sum(log(1+exp(Xcat%*%c(beta)))*counts))/sum(counts)
        #   } 
        #   fit.AD <- optim(par = rep(0, ncol(Xcat)), logLik_AD, method = "BFGS") 
        #   se <- sqrt(diag(solve(hessian(func = function(x) logLik_AD(x)*sum(counts), x = fit.AD$par))))
        #   res <- data.frame(est = fit.AD$par, se = se)
        #   rownames(res) <- colnames(Xcat) 
        #   control$final_output = res 
        # }else if(control$model == "OLGLMM"){
        #   for(site_i in control$site[1]){
        #     AD = pdaGet(paste0(site_i,'_initialize'),config)
        #     Xtable <- AD$Xtable
        #     SY <- Xtable$SY
        #     count <- Xtable$n
        #     new_Xtable <- Xtable[,1:(1+length(control$variables))] 
        #     new_X <- new_Xtable[rep(seq_len(nrow(new_Xtable)), times = count), ]
        #     generate_Y <- function(Y_count, overall_count){
        #       sub_Y <- rep(c(1,0), times = c(Y_count, overall_count- Y_count))
        #     }
        #     new_Y <- unlist(mapply(generate_Y,SY,overall_count = count)) 
        #     output_0 = cbind(new_Y, new_X)
        #     colnames(output_0) = c(control$outcome,"intercept",control$variables)
        #     output_0$site <- site_i
        #   }
        #   
        #   for(site_i in control$sites[-1]){
        #     AD = pdaGet(paste0(site_i,'_initialize'),config)
        #     Xtable <- AD$Xtable
        #     SY <- Xtable$SY
        #     count <- Xtable$n
        #     new_Xtable <- Xtable[,1:(1+length(control$variables))] 
        #     new_X <- new_Xtable[rep(seq_len(nrow(new_Xtable)), times = count), ]
        #     generate_Y <- function(Y_count, overall_count){
        #       sub_Y <- rep(c(1,0), times = c(Y_count, overall_count- Y_count))
        #     }
        #     new_Y <- unlist(mapply(generate_Y,SY,overall_count = count)) 
        #     output = cbind(new_Y, new_X)
        #     colnames(output) = c(control$outcome,"intercept",control$variables)
        #     output$site <- site_i
        #     output_0 <- rbind(output_0, output)
        #   } 
        #   
        #   fit.pool.recons <- MASS::glmmPQL(as.formula(paste(control$outcome, paste(control$variables, collapse = "+"), sep = '~')), 
        #                              ~1|site, 
        #                              data = output_0,
        #                              family='binomial')
        #   control$est <- fit.pool.recons$coefficients
        #   control$varFix <- fit.pool.recons$varFix
        #   control$rand_se <- fit.pool.recons$sigma 
      }else if(control$model=='COLA'){
        # COLA does not need derivative, thus estimate after initialize...
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
        }else{ # all other ODAX..
          init_i = pdaGet(paste0(control$sites[1],'_initialize'),config)
          bhat <-init_i$bhat_i 
          vbhat <- init_i$Vhat_i
          for(site_i in control$sites[-1]){
            init_i <- pdaGet(paste0(site_i,'_initialize'),config)
            bhat = rbind(bhat, init_i$bhat_i)
            vbhat = rbind(vbhat, init_i$Vhat_i)
          }
        }
        
        site_size = c()
        for(site_i in control$sites){ 
          init_i <- pdaGet(paste0(site_i,'_initialize'),config)
          site_size <- c(site_size, init_i$site_size)
        }
        # print(site_size)
        
        ## estimate for pda init: meta, or median, or lead est?...
        if(control$init_method == 'meta'){
          binit = apply(as.data.frame(bhat/vbhat),2,function(x){sum(x, na.rm = TRUE)})/apply(as.data.frame(1/vbhat),2,function(x){sum(x, na.rm = TRUE)})
          # vinit = 1/apply(1/vbhat,2,function(x){sum(x, na.rm = TRUE)}) 
          mymessage('meta (inv var weighted avg) as initial est:')
        } else if(control$init_method == 'median'){ 
          binit = apply(bhat, 2, median, na.rm=T) 
          mymessage('median as initial est:')
        # } else if(control$init_method == 'weighted.median'){
        #   binit = apply(bhat, 2, function(x) weighted.median(x, site_size))
        #   mymessage('median (site size weighted) as initial est:')
        } else if(control$init_method == 'lead'){
          binit = bhat[control$sites==control$lead_site,]
          mymessage('lead site est as initial est:')
        } #print(res)
        
        control$beta_init = binit
        if(any(is.na(control$beta_init))){
          control$beta_init[is.na(control$beta_init)] = 0
          warning('some coefs are all NAs and use 0 as init est, use caution and check your X variables!!!')
        }
        
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
        mymessage("You are done!")
      } 
      
      ## robust var est for ODACH_CC est (might also for other ODAX later?)
      if (control$model == "ODACH_CC"){
        px = length(control$beta_init) 
        K = length(control$sites)
        S_i_sum = array(NA,c(px,px, K))
        for(i in 1:K){ 
          site_i = control$sites[i]
          init_i <- pdaGet(paste0(site_i,'_derive'),config)
          S_i_sum[,,i] <- init_i$S_i
        }
        control$S_i_sum = apply(S_i_sum,c(1,2),sum,na.rm=T) 
        # print(control$S_i_sum)
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
  
  ## update control with next step
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
  
  mymessage(mes)
  pdaPut(control,'control',config,upload_without_confirm,silent_message,digits)
  
  control
}
