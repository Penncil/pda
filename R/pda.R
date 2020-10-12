# Copyright 2020 Penn Computing Inference Learning (PennCIL) lab
#
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
#' @usage pdaPut(obj,name,config)
#' @param obj R object to encode as json and uploaded to cloud
#' @param name of file
#' @param config a list of variables for cloud configuration
#' @importFrom utils menu
#' @return NONE
#' @seealso \code{pda}
#' @export
pdaPut <- function(obj,name,config){
    obj_Json <- jsonlite::toJSON(obj)
    file_name <- paste0(name, '.json')
    if(interactive()) {
      print(paste("Put",file_name,"on public cloud:"))
      print(obj_Json)
      authorize = menu(c("Yes", "No"), title="Allow upload?")
    } else {
      authorize = "1"
    }
    if (authorize != 1) {
      print("file not uploaded.")
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
        print(paste("putting:",url))
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
  if(!is.null(dir)) {
      config$dir = dir
  } else if (pda_dir!='') {
      config$dir = pda_dir
  }
  if(!is.null(uri)) {
      config$uri = uri
  } else if (pda_uri!='') {
      config$uri = pda_uri
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
#' @return control
#' @seealso \code{pdaPut}, \code{pdaList}, \code{pdaGet}, \code{getCloudConfig} and \code{pdaSync}.
#' @import stats survival rvest jsonlite data.table httr Rcpp RcppArmadillo
#'          
#' @references
#' Michael I. Jordan, Jason D. Lee & Yun Yang (2019) Communication-Efficient Distributed Statistical Inference, \cr
#'  \emph{Journal of the American Statistical Association}, 114:526, 668-681 \cr 
#'  \url{https://doi.org/10.1080/01621459.2018.1429274}.\cr 
#' (ODAL) Rui Duan, et al. (2020) Learning from electronic health records across multiple sites: \cr 
#'  A communication-efficient and privacy-preserving distributed algorithm. \cr 
#'  \emph{Journal of the American Medical Informatics Association}, 27.3:376–385,
#'  \cr \url{https://doi.org/10.1093/jamia/ocz199}.\cr 
#' (ODAC) Rui Duan, et al. (2020) Learning from local to global: An efficient distributed algorithm for modeling time-to-event data. \cr
#'   \emph{Journal of the American Medical Informatics Association}, 27.7:1028–1036, \cr 
#'    \url{https://doi.org/10.1093/jamia/ocaa044}.
#' @examples
#' require(survival)
#' require(data.table)
#' require(pda)
#' data(lung)
#' 
#' # Create a number of sites, split the lung data amongst them
#' sites = c('site1', 'site2', 'site3')
#' set.seed(42)
#' lung2 <- lung[,2:5]
#' lung2$sex <- lung2$sex - 1
#' lung2$status <- ifelse(lung2$status == 2, 1, 0)
#' lung_split <- split(lung2, sample(1:length(sites), nrow(lung), replace=TRUE))
#' ## fit logistic reg using pooled data
#' fit.pool <- glm(status ~ age + sex, family = 'binomial', data = lung2)
#' 
#' ## In the example below we aim to use PDA ODAL to obtain a surrogate estimator that is 
#' ## close to the pooled estimate. Accounts (site1, site2, site3) and password 
#' ## (WLjySoaZnSqMNswowg) are given to the 3 example sites at the server https://pda.one. 
#' ## Each site can access via web browser to check the communication of the summary stats.
#' 
#' # ############################  STEP 1: initialize  ###############################
#' ## lead site1: please review and enter "1" to allow putting the control file to the server
#' control <- list(project_name = 'Lung cancer study',
#'                 step = 'initialize',     
#'                 sites = sites,
#'                 heterogeneity = FALSE,
#'                 model = 'ODAL',
#'                 outcome = "status",
#'                 variables = c('age', 'sex'),
#'                 optim_maxit = 100,
#'                 lead_site = sites[1],
#'                 upload_date = as.character(Sys.time()) )
#' Sys.setenv(PDA_USER = 'site1', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
#' pda(site_id = 'site1', control = control)  
#' # now the group would see control.json	at https://pda.one
#' 
#' # remote site3: please review and enter "1" to allow putting your local estimate to the server  
#' i <- 3
#' Sys.setenv(PDA_USER = 'site3', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
#' control <- pda(ipdata = lung_split[[i]], site_id = sites[i])
#' # now the group would see site3_initialize.json	at https://pda.one
#' 
#' ## remote site2: please review and enter "1" to allow putting your local estimate to the server
#' i <- 2
#' Sys.setenv(PDA_USER = 'site2', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
#' control <- pda(ipdata = lung_split[[i]], site_id = sites[i])
#' # now the group would see site2_initialize.json	at https://pda.one
#' 
#' ## lead site1: please review and enter "1" to allow putting your local estimate to the server
#' i <- 1
#' Sys.setenv(PDA_USER = 'site1', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
#' control <- pda(ipdata = lung_split[[i]], site_id = sites[i]) 
#' # now the group would see site1_initialize.json	at https://pda.one
#' # control.json is also automatically updated  
#' 
#' ## if lead site1 initialized before other sites, 
#' ## lead site1: uncomment to synchoronize the control before STEP 2
#' # Sys.setenv(PDA_USER = 'site1', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
#' # pda(site_id = 'site1', control = control)
#' # config <- getCloudConfig(site_id = 'site1')
#' # pdaSync(config)
#' 
#' # ############################  STEP 2: derivative  ###############################
#' ## remote site3: please review and enter "1" to allow putting your derivatives to the server   
#' i <- 3
#' Sys.setenv(PDA_USER = 'site3', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
#' control <- pda(ipdata = lung_split[[i]], site_id = sites[i])
#' # now the group would see site3_derive.json	at https://pda.one
#' 
#' ## remote site2: please review and enter "1" to allow putting your derivatives to the server    
#' i <- 2
#' Sys.setenv(PDA_USER = 'site2', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
#' control <- pda(ipdata = lung_split[[i]], site_id = sites[i])
#' # now the group would see site2_derive.json	at https://pda.one
#' 
#' ## lead site1: please review and enter "1" to allow putting your derivatives to the server      
#' i <- 1
#' Sys.setenv(PDA_USER = 'site1', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
#' control <- pda(ipdata = lung_split[[i]], site_id = sites[i])
#' # now the group would see site1_derive.json	at https://pda.one
#' 
#' # ############################  STEP 3: estimate  ###############################
#' ## lead site1: 
#' ## please review and enter "1" to allow putting the surrogate estimate to the server     
#' i <- 1
#' Sys.setenv(PDA_USER = 'site1', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
#' control <- pda(ipdata = lung_split[[i]], site_id = sites[i])
#' # now the group would see site1_estimate.json	at https://pda.one
#' ## compare the surrogate estimate with the pooled estimate
#' config <- getCloudConfig(site_id = 'site1')
#' fit.odal <- pdaGet(name = 'site1_estimate', config = config)
#' cbind(b.pool=fit.pool$coef, 
#'       b.odal=fit.odal$btilde, 
#'       sd.pool=summary(fit.pool)$coef[,2], 
#'       sd.odal=sqrt(diag(solve(fit.odal$Htilde)/nrow(lung2))))
#' ## the PDA ODAL is now completed! 
#' ## All the sites can still run their surrogate estimates and broadcast them. 
#' 
#' ## remote site2: (optional)  
#' i <- 2
#' Sys.setenv(PDA_USER = 'site2', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
#' control <- pda(ipdata = lung_split[[i]], site_id = sites[i])
#' 
#' ## remote site3: (optional)
#' i <- 3
#' Sys.setenv(PDA_USER = 'site3', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
#' control <- pda(ipdata = lung_split[[i]], site_id = sites[i]) 
#' 
#' 
#' ## If all the sites broadcast their surrogate estimates, 
#' ## a final synthesize step can further improve the estimate.
#' ## lead site1: uncomment to synchoronize the control before STEP 4
#' Sys.setenv(PDA_USER = 'site1', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
#' pda(site_id = 'site1', control = control)
#' config <- getCloudConfig(site_id = 'site1')
#' pdaSync(config)
#' 
#' # ########################  STEP 4: synthesize (optional)  ######################## 
#' ## lead site1:     
#' i <- 1
#' Sys.setenv(PDA_USER = 'site1', PDA_SECRET = 'WLjySoaZnSqMNswowg', PDA_URI = 'https://pda.one')
#' control <- pda(ipdata = lung_split[[i]], site_id = sites[i])
#'
#' @return control
#' @export
pda <- function(ipdata=NULL,site_id,control=NULL,dir=NULL,uri=NULL,secret=NULL){
                config<-getCloudConfig(site_id,dir,uri,secret)
  #add a control if one was provided
  if(!(is.null(control)) &&  config$site_id==control$lead_site) { # control$sites[1]
           pdaPut(obj=control,name='control',config=config)
           return(control)    # ?
  }
  control = pdaGet('control',config)
  cat('You are performing Privacy-preserving Distributed Algorithm (PDA, https://github.com/Penncil/pda): \n')
  cat('your site = ', config$site_id, '\n')
  n = nrow(ipdata)
  formula<-as.formula(
    paste(control$outcome,
    paste(control$variables, collapse = " + "),
  sep = ' ~'))
  mf = model.frame(formula, ipdata)
  
  if(control$model=='ODAL'){
    ODAL.steps<-c('initialize','derive','estimate','synthesize')
    ODAL.family<-'binomial'
  }else if(control$model=='ODAC'){
    ODAC.steps<-c('initialize','derive','derive_UWZ','estimate','synthesize')
    ODAC.family<-'cox'
  }
  
  family = get(paste0(control$model,'.family'))
  if(family=='cox'){  
    ipdata = data.table::data.table(time=as.numeric(model.response(mf))[1:n], 
                        status=as.numeric(model.response(mf))[-c(1:n)], 
                        model.matrix(formula, mf)[,-1])
    control$risk_factor = colnames(ipdata)[-c(1:2)]
  }else{
    ipdata = data.table::data.table(status=as.numeric(model.response(mf)), 
                        model.matrix(formula, mf))
    control$risk_factor = colnames(ipdata)[-1]
  }
  if(is.character(control$step)){
    step_function<-paste0(control$model,'.',control$step)
    step_obj<-get(step_function)(ipdata, control, config)
    if(control$step=='estimate'){
      print("Congratulations, the PDA is completed! You can continue broadcasting your surrogate estimate to further synthesize them.")
    }
    pdaPut(step_obj,paste0(config$site_id,'_',control$step),config)
    #sync needed?
    if(config$site_id==control$lead_site) {
           control<-pdaSync(config)
    }
  }
  control
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
pdaSync <- function(config){  
  control = pdaGet('control',config)
  if(control$model=='ODAL'){
    ODAL.steps<-c('initialize','derive','estimate','synthesize')
    ODAL.family<-'binomial'
  }else if(control$model=='ODAC'){
    ODAC.steps<-c('initialize','derive','derive_UWZ','estimate','synthesize')
    ODAC.family<-'cox'
  }
  
  files<-pdaList(config) 
  if(all(paste0(control$sites,"_",control$step) %in% files)){
    if(control$step=="initialize"){
      init_i <- pdaGet(paste0(control$lead_site,'_initialize'),config)
      bhat <-init_i$bhat_i 
      vbhat <- init_i$Vhat_i
      for(site_i in control$sites){
        if(site_i!=control$lead_site){
          init_i <- pdaGet(paste0(site_i,'_initialize'),config)
          bhat = rbind(bhat, init_i$bhat_i)
          vbhat = rbind(vbhat, init_i$Vhat_i)
        }
      }
      #estimate from meta-analysis
      bmeta = apply(bhat/vbhat,2,function(x){sum(x, na.rm = TRUE)})/apply(1/vbhat,2,function(x){sum(x, na.rm = TRUE)})
      vmeta = 1/apply(1/vbhat,2,function(x){sum(x, na.rm = TRUE)})
      res = list(bmeta = bmeta, vmeta = vmeta)
      cat('meta analysis (inverse variance weighted average) result:')
      #print(res)
      control$beta_init = bmeta
      mes <- 'beta_init added, step=2 (derivatives)! \n'
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
    cat(mes)
    pdaPut(control,'control',config)
    control
  } 
}
