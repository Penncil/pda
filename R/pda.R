# https://style.tidyverse.org/functions.html#naming

# require(survival)
# require(data.table)
# Rcpp::sourceCpp('src/rcpp_coxph.cpp')
## broadcast: upload/download shared info to/from the cloud folder

#' @useDynLib pda
#' @title Function to upload object to cloud as json
#' 
#' @usage pdaPut(obj,name,cloud_config)
#' @author Chongliang Luo, Steven Vitale
#' @param obj R object to encode as json and uploaded to cloud
#' @param name of file
#' @param cloud_config a list of variables for cloud configuration
#' @export
#' @return  
pdaPut <- function(obj,name,cloud_config){
    obj_Json <- jsonlite::toJSON(obj)
    file_name <- paste0(name, '.json')
    #print(paste("Put",file_name,"on public cloud:"))
    #print(obj_Json)
    if(interactive()) {
      authorize = menu(c("Yes", "No"), title="Allow upload?")
    } else {
      authorize = "1"
    }
    if (authorize != 1) {
      print("file not uploaded.")
      return(FALSE)
    }
    # the file to upload
    if (is.character(cloud_config$dir)) {
        file_path <- paste0(cloud_config$dir,'/', file_name)
    } else {
        file_path <- paste0(tempdir(),'/', file_name)
    }
    #print(paste0("writing",name,"to",file_path))
    write(obj_Json, file_path)
    if (is.character(cloud_config$uri)) {
        # create the url target of the file
        url <- file.path(cloud_config$uri, file_name)
        # webdav PUT request to send a file to cloud
        r<-httr::PUT(url, body = httr::upload_file(file_path), httr::authenticate(cloud_config$site_id, cloud_config$secret, 'digest'))
        print(paste("putting:",url))
    }
}


#' @useDynLib pda
#' @title Function to download json and return as object
#' @usage pdaGet(name)
#' @author Chongliang Luo, Steven Vitale
#' @param name of file
#' @return  
#' @export
pdaGet <- function(name,cloud_config){
    file_name <- paste0(name, '.json')
    #print(paste("Get",file_name,"from public cloud:"))
    # the file to upload
    if (is.character(cloud_config$dir)) {
        file_path <- paste0(cloud_config$dir,'/', file_name)
    } else {
        file_path <- paste0(tempdir(),'/', file_name)
    }
    if (is.character(cloud_config$uri)) {
        url <- file.path(cloud_config$uri, file_name)
        #write the file from GET request to file_path
        res<-httr::GET(url, httr::write_disk(file_path, overwrite = TRUE), httr::authenticate(cloud_config$site_id, cloud_config$secret, 'digest'))
        #print(paste("getting:",url))
    } 
    obj<-jsonlite::fromJSON(file_path)
    return(obj)
}


#' @useDynLib pda
#' @title gather cloud settings into one list
#' @usage getCloudConfig()
#' @author Chongliang Luo, Steven Vitale
#' @param name of file
#' @export
getCloudConfig <- function(site_id=NULL){
  cloud_config=list()
  pda_user<-Sys.getenv('PDA_USER')
  pda_secret<-Sys.getenv('PDA_SECRET')
  pda_uri<-Sys.getenv('PDA_URI')
  pda_dir<-Sys.getenv('PDA_DIR')
  if(!is.null(site_id)) {
    cloud_config$site_id=site_id
  } else if(pda_user!='') {
    cloud_config$site_id=pda_user
  }
  if (pda_secret!='') cloud_config$secret = pda_secret
  if (pda_uri!='') cloud_config$uri = pda_uri
  if (pda_dir!='') cloud_config$dir = pda_dir
  return(cloud_config);
}

#' @useDynLib pda
#' @title Function to download json and return as object
#' @param control pda control object
#' @param site_id name of site
#' @export
pdaInit <- function(control,site_id=NULL){
  cloud_config<-getCloudConfig(site_id)
  pdaPut(control,'control',cloud_config)
}

#' @useDynLib pda
#' @import stats
#' @import Rcpp
#' @import survival 
#' @title PDA: Privacy-preserving Distributed Algorithm
#' 
#' @description  Fit Privacy-preserving Distributed Algorithms for linear, logistic, 
#'                Poisson and Cox PH regression with possible heterogeneous data across sites.
#' @usage pda(data = ipdata, mysite = NULL)
#' @author Chongliang Luo, Steven Vitale, Jiayi Tong, Rui Duan, Yong Chen
#' 
#' @param ipdata  Local IPD data in data frame, should include at least one column for the outcome and one column for the covariates 
#' @param mysite Character, site name as decided by all the sites when initialize the collaboration
#'
#'          
#' @references Jordan, Michael I., Jason D. Lee, and Yun Yang. "Communication-efficient distributed statistical inference." JASA (2018).
#' @examples
#'  ## Steve can you make an ODAl example here as we tested?
#' 
#' @export
pdaStep <- function(ipdata,site_id=NULL){
  cloud_config<-getCloudConfig(site_id)
  control = pdaGet('control',cloud_config)
  cat('You are performing Privacy-preserving Distributed Algorithm (PDA, https://github.com/Penncil/pda): \n')
  cat('your site = ', cloud_config$site_id, '\n')
  n = nrow(ipdata)
  formula<-as.formula(
    paste(control$outcome,
    paste(control$variables, collapse = " + "),
  sep = ' ~'))
  mf = model.frame(formula, ipdata)
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
    # # family = 
    # # gaussian: DLM / DLMM, 
    # # binomial: ODAL / ODALH / ODALR, 
    # # poisson:  ODAP / ODAH / DPLR, 
    # # cox:      ODAC / ODACH
    # if(control$heterogeneity==FALSE) control$model = 'ODAL'
    # else control$model = 'ODALH'
  }
  if(is.character(control$step)){
    step_function<-paste0(control$model,'.',control$step)
    step_obj<-get(step_function)(ipdata, control, cloud_config)
    pdaPut(step_obj,paste0(cloud_config$site_id,'_',control$step),cloud_config)
  } else {
   print('finished')
  }
}


#' update the PDA control, used by the master site
#' @usage pdaSync <- function(config,site_id)
#' 
#' @param control_update Logical, update the PDA control?
#' 
#' @export
pdaSync <- function(cloud_config,site_id=NULL){
  cloud_config<-getCloudConfig(site_id)
  control = pdaGet('control',cloud_config)
  if(control$step=="initialize"){
    init_i <- pdaGet(paste0(control$lead_site,'_initialize'),cloud_config)
    bhat <-init_i$bhat_i 
    vbhat <- init_i$Vhat_i
    for(site_i in control$sites){
      if(site_i!=control$lead_site){
        init_i <- pdaGet(paste0(site_i,'_initialize'),cloud_config)
        bhat = rbind(bhat, init_i$bhat_i)
        vbhat = rbind(vbhat, init_i$Vhat_i)
      }
    }
   
    #estimate from meta-analysis
    bmeta = apply(bhat/vbhat,2,function(x){sum(x, na.rm = T)})/apply(1/vbhat,2,function(x){sum(x, na.rm = T)})
    vmeta = 1/apply(1/vbhat,2,function(x){sum(x, na.rm = T)})
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
  pdaPut(control,'control',cloud_config)
  return(control)
}
