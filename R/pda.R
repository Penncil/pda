#' @useDynLib pda
#' @title Function to upload object to cloud as json
#' 
#' @usage pdaPut(obj,name,config)
#' @param obj R object to encode as json and uploaded to cloud
#' @param name of file
#' @param config a list of variables for cloud configuration
#' @importFrom utils menu
#' @return NONE
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
#' @return list of files
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
#' @return data from json file
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
#' @return list of cloud parameters
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
#' @import stats survival rvest jsonlite data.table httr Rcpp RcppArmadillo
#'
#'          
#' @references Jordan, Michael I., Jason D. Lee, and Yun Yang. "Communication-efficient distributed statistical inference." JASA (2018). https://doi.org/10.1080/01621459.2018.1429274
#'             Rui Duan, et al. "Learning from electronic health records across multiple sites: A communication-efficient and privacy-preserving distributed algorithm". Journal of the American Medical Informatics Association, 2020, https://doi.org/10.1093/jamia/ocz199
#'             Rui Duan, et al. "Learning from local to global: An efficient distributed algorithm for modeling time-to-event data". Journal of the American Medical Informatics Association, 2020, https://doi.org/10.1093/jamia/ocaa044
#' @examples
#'require(survival)
#'require(data.table)
#'require(pda)
#'data(lung)
#'#create a number of sites, split the lung data amongst them
#' sites = c('site1', 'site2', 'site3')
#' set.seed(42)
#' lung2<-lung[,2:5]
#' lung2$sex <- lung2$sex-1
#' lung2$status <- ifelse(lung2$status == 2, 1, 0)
#' lung_split<-split(lung2, sample(1:length(sites), nrow(lung), replace=TRUE))
#' 
#' control <- list(project_name = 'Lung cancer study',
#'                 step = 'initialize',    #' current step, updated by lead
#'                 sites = sites,
#'                 heterogeneity = FALSE,
#'                 model = 'ODAL',
#'                 outcome = "status",
#'                 variables = c('age', 'sex'),
#'                 optim_maxit=100,
#'                 lead_site = sites[1],
#'                 upload_date = as.character(Sys.time()) )
#' 
#' #RUN BY LEAD ONLY , check results at https://pda.one/003   # , PDA_DIR='/Users/chl18019/Dropbox/PDA-git'
#' Sys.setenv(PDA_USER='site1', PDA_SECRET='WLjySoaZnSqMNswowg', PDA_URI='https://pda.one/003')
#' pda(site_id='site1',control=control)
#' 
#' # config<-getCloudConfig(site_id='site1' )
#' 
#' #run pda until step is empty
#' while (is.character(control$step)) {
#'   print(paste("step:",control$step))
#'   #cycle through sites
#'   for(i in length(sites):1) {      # , PDA_DIR='/Users/chl18019/Dropbox/PDA-git'
#'     Sys.setenv(PDA_USER=paste0('site',i), PDA_SECRET='WLjySoaZnSqMNswowg', PDA_URI='https://pda.one/003')
#'     control<-pda(ipdata=lung_split[[i]],site_id=sites[i])
#'   }
#' }
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
    pdaPut(step_obj,paste0(config$site_id,'_',control$step),config)
    #sync needed?
    if(config$site_id==control$sites[1]) {
           control<-pdaSync(config)
    }
  }
  control
}


#' @useDynLib pda
#' @title pda synchronize 
#' 
#' @description  update pda control if ready (run by lead)
#' @usage pdaSync(config)
#' @param config cloud configuration
#' @return control
#'   
## maybe avoid this step to give lead site more freedom to manually adjust the procedure 
# for example, exclude sites that are not valid after initialization, manually specify a init value
pdaSync <- function(config){  
  control = pdaGet('control',config)
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
