# Copyright 2020 Penn Computing Inference Learning (PennCIL) lab
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
      if(!is.null(config$uri)){
        print(paste("Put",file_name,"on public cloud:"))
      }else{
        print(paste("Put",file_name,"on local directory", config$dir, ':'))
      }
      print(obj_Json)
      authorize = menu(c("Yes", "No"), title="Allow transfer?")
    } else {
      authorize = "1"
    }
    if (authorize != 1) {
      print("file not transferred. You can specify file transfer setting in pda().")
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

  if(!is.null(uri)) {
      config$uri = uri
  } else if (pda_uri!='') {
      config$uri = pda_uri
  } else{
    print('no cloud uri found! ')
  }
  
  if(!is.null(dir)) {
    config$dir = dir
  } else if (pda_dir!='') {
    config$dir = pda_dir
  }else{
    cat('no public or local directory supplied, use local temporary:', tempdir())
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
#' ## In the toy example below we aim to analyze the association of lung status with 
#' ## age and sex using logistic regression, data(lung) from 'survival', we randomly 
#' ## assign to 3 sites: 'site1', 'site2', 'site3'. we demonstrate using PDA ODAL can 
#' ## obtain a surrogate estimator that is close to the pooled estimate. We run the 
#' ## example in local directory. In actual collaboration, account/password for pda server 
#' ## will be assigned to the sites at the server https://pda.one.
#' ## Each site can access via web browser to check the communication of the summary stats.
#' 
#' ## for more examples, see demo(ODAC_lung_cancer) and demo(ODAP_CrabSatellites)
#' 
#' # Create 3 sites, split the lung data amongst them
#' sites = c('site1', 'site2', 'site3')
#' set.seed(42)
#' lung2 <- lung[,c('time', 'status', 'age', 'sex')]
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
#' ## assume lead site1: enter "1" to allow transferring the control file  
#' pda(site_id = 'site1', control = control, dir = getwd())
#' ## in actual collaboration, account/password for pda server will be assigned, thus:
#' # pda(site_id = 'site1', control = control, uri = 'https://pda.one', secret='abc123')
#' ## you can also set your environment variables, and no need to specify them in pda:
#' # Sys.setenv(PDA_USER = 'site1', PDA_SECRET = 'abc123', PDA_URI = 'https://pda.one')
#' # pda(site_id = 'site1', control = control)
#' 
#' ##' assume remote site3: enter "1" to allow tranferring your local estimate 
#' pda(site_id = 'site3', ipdata = lung_split[[3]], dir=getwd())
#' 
#' ##' assume remote site2: enter "1" to allow tranferring your local estimate  
#' pda(site_id = 'site2', ipdata = lung_split[[2]], dir=getwd())
#' 
#' 
#' ##' assume lead site1: enter "1" to allow tranferring your local estimate  
#' ##' control.json is also automatically updated
#' pda(site_id = 'site1', ipdata = lung_split[[1]], dir=getwd())
#' 
#' ##' if lead site1 initialized before other sites,
#' ##' lead site1: uncomment to sync the control before STEP 2
#' #' pda(site_id = 'site1', control = control)
#' #' config <- getCloudConfig(site_id = 'site1')
#' #' pdaSync(config)
#' 
#' #' ############################'  STEP 2: derivative  ############################ 
#' ##' assume remote site3: enter "1" to allow tranferring your derivatives  
#' pda(site_id = 'site3', ipdata = lung_split[[3]], dir=getwd())
#' 
#' ##' assume remote site2: enter "1" to allow tranferring your derivatives  
#' pda(site_id = 'site2', ipdata = lung_split[[2]], dir=getwd())
#' 
#' ##' assume lead site1: enter "1" to allow tranferring your derivatives  
#' pda(site_id = 'site1', ipdata = lung_split[[1]], dir=getwd())
#' 
#' 
#' #' ############################'  STEP 3: estimate  ############################ 
#' ##' assume lead site1: enter "1" to allow tranferring the surrogate estimate  
#' pda(site_id = 'site1', ipdata = lung_split[[1]], dir=getwd())
#' 
#' ##' the PDA ODAL is now completed!
#' ##' All the sites can still run their own surrogate estimates and broadcast them.
#' 
#' ##' compare the surrogate estimate with the pooled estimate
#' config <- getCloudConfig(site_id = 'site1', dir=getwd())
#' fit.odal <- pdaGet(name = 'site1_estimate', config = config)
#' cbind(b.pool=fit.pool$coef,
#'       b.odal=fit.odal$btilde,
#'       sd.pool=summary(fit.pool)$coef[,2],
#'       sd.odal=sqrt(diag(solve(fit.odal$Htilde)/nrow(lung2))))
#' 
#' @return control
#' @export
pda <- function(ipdata=NULL,site_id,control=NULL,dir=NULL,uri=NULL,secret=NULL){
  config <- getCloudConfig(site_id,dir,uri,secret)
  #add a control if one was provided
  if(!(is.null(control)) &&  config$site_id==control$lead_site) { # control$sites[1]
           pdaPut(obj=control,name='control',config=config)
           return(control)    # ?
  }
  control = pdaGet('control',config)
  cat('You are performing Privacy-preserving Distributed Algorithm (PDA, https://github.com/Penncil/pda): \n')
  cat('your site = ', config$site_id, '\n')
 
  if(control$model=='ODAL'){
    ODAL.steps<-c('initialize','derive','estimate','synthesize')
    ODAL.family<-'binomial'
  }else if(control$model=='ODAP'){
    ODAP.steps<-c('initialize','derive','estimate','synthesize')
    # for ODAP need to specify family (poisson, ztpoisson, quasipoisson, ztquasipoisson, hurdle) in control
    ODAP.family<-control$family  
  }else if(control$model=='ODAC'){
    ODAC.steps<-c('initialize','derive','derive_UWZ','estimate','synthesize')
    ODAC.family<-'cox'
  }

  family = get(paste0(control$model,'.family'))
  n = nrow(ipdata)
  if(family=='hurdle'){           # count and zero parts for hurdle, Xcount first
    variables <- control$variables_hurdle_count
  }else{
    variables <- control$variables
  }
  formula <- as.formula(paste(control$outcome, paste(variables, collapse = "+"), sep = '~'))
  mf <- model.frame(formula, ipdata)
  
  # create ipdata via model.matrix to make dummy variables for categorical covariates...
  if(control$model=='ODAC'){  
    ipdata = data.table::data.table(time=as.numeric(model.response(mf))[1:n], 
                        status=as.numeric(model.response(mf))[-c(1:n)], 
                        model.matrix(formula, mf)[,-1])
    control$risk_factor = colnames(ipdata)[-c(1:2)]
  }else if(control$model=='ODAL'){
    ipdata = data.table::data.table(status=as.numeric(model.response(mf)), 
                                    model.matrix(formula, mf))
    control$risk_factor = colnames(ipdata)[-1]
  }else if(control$model=='ODAP'){
    if(family=='hurdle'){        # count and zero parts for hurdle
      X_count = data.table::data.table(model.matrix(formula, mf))
      # also make design X_zero
      formula <- as.formula(paste(control$outcome, paste(control$variables_hurdle_zero, collapse = "+"), sep = '~'))
      mf <- model.frame(formula, ipdata)
      X_zero = data.table::data.table(model.matrix(formula, mf))
      ipdata <- list(ipdata=ipdata, X_count=X_count, X_zero=X_zero, 
                     offset = ifelse(is.character(control$offset), ipdata[,control$offset], 0))  
      control$risk_factor = c('Intercept', control$variables_hurdle_count, 'Intercept', control$variables_hurdle_zero)   
      # colnames(ipdata)[-c(1:2)]
    }else{
      ipdata = data.table::data.table(outcome=as.numeric(model.response(mf)), 
                                      offset=ifelse(is.character(control$offset), ipdata[,control$offset], 0),
                                      model.matrix(formula, mf))
      control$risk_factor = colnames(ipdata)[-c(1:2)]
    }
  }
  
  if(is.character(control$step)){
    step_function <- paste0(control$model,'.',control$step)
    step_obj <- get(step_function)(ipdata, control, config)
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
  }else if(control$model=='ODAP'){
    ODAP.steps<-c('initialize','derive','estimate','synthesize')
    # for ODAP need to specify family (poisson, ztpoisson, quasipoisson, ztquasipoisson, hurdle) in control
    ODAP.family<-control$family  
  }else if(control$model=='ODAC'){
    ODAC.steps<-c('initialize','derive','derive_UWZ','estimate','synthesize')
    ODAC.family<-'cox'
  }
  
  files<-pdaList(config) 
  if(all(paste0(control$sites,"_",control$step) %in% files)){
    if(control$step=="initialize"){
      init_i <- pdaGet(paste0(control$lead_site,'_initialize'),config)
      if(control$family=='hurdle'){
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
        cat('meta analysis (inverse variance weighted average) result:')
        #print(res)
        control$beta_zero_init = bmeta_zero
        control$beta_count_init = bmeta_count
      }else{
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
      }
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
