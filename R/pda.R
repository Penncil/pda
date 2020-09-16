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
#' @importFrom utils menu
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
#' @title Function to list available objects
#' @usage pdaList(name)
#' @author Chongliang Luo, Steven Vitale
#' @importFrom rvest html_nodes
#' @import httr 
#' @return  
pdaList <- function(cloud_config){
    if (is.character(cloud_config$uri)) {
        res<-httr::GET(cloud_config$uri, httr::authenticate(cloud_config$site_id, cloud_config$secret, 'digest'))
        files<-rvest::html_nodes(httr::content(res), xpath = "//table//td/a") 
        files<-files[lapply(files,length)>0]
        files<-regmatches(files, gregexpr("(?<=\")(.*?)(json)(?=\")", files, perl = TRUE))
    } else if (is.character(cloud_config$dir)) {
        files<-list.files(cloud_config$dir,pattern = "\\.json$") 
    } else {
        files<-list.files(tempdir(),pattern = "\\.json$") 
    }
    files<-substr(files,1,nchar(files)-5)
    return(files)
}


#' @useDynLib pda
#' @title Function to download json and return as object
#' @usage pdaGet(name)
#' @author Chongliang Luo, Steven Vitale
#' @param name of file
#' @return  
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
#' @usage getCloudConfig()
#' @author Chongliang Luo, Steven Vitale
#' @param site_id site identifier
#' @param dir shared directory path if using flat files
#' @param uri web uri if using web service
#' @param secret web token if using web service
getCloudConfig <- function(site_id,dir=NULL,uri=NULL,secret=NULL){
  cloud_config<-list()
  cloud_config$site_id=site_id
  if(!is.null(secret)) {
      cloud_config$secret = secret
  }
  if(!is.null(dir)) {
      cloud_config$dir = dir
  }
  if(!is.null(uri)) {
      cloud_config$uri = uri
  }
  cloud_config;
}

#' @useDynLib pda
#' @title PDA: Privacy-preserving Distributed Algorithm
#' 
#' @description  Fit Privacy-preserving Distributed Algorithms for linear, logistic, 
#'                Poisson and Cox PH regression with possible heterogeneous data across sites.
#' @usage pda(data = ipdata, mysite = NULL)
#' @author Chongliang Luo, Steven Vitale, Jiayi Tong, Rui Duan, Yong Chen
#' 
#' @param ipdata  Local IPD data in data frame, should include at least one column for the outcome and one column for the covariates 
#' @param site_id Character site name
#' @param control pda control data
#'
#'          
#' @references Jordan, Michael I., Jason D. Lee, and Yun Yang. "Communication-efficient distributed statistical inference." JASA (2018).
#' @examples
#'require(survival)
#'require(data.table)
#'require(pda)
#'data(lung)
#'#create a number of sites, split the lung data amongst them
#'sites = c('site1', 'site2', 'site3')
#'set.seed(42)
#'lung2<-lung[,2:5]
#'lung2$sex <- lung2$sex-1
#'lung2$status <- ifelse(lung2$status == 2, 1, 0)
#'lung_split<-split(lung2, sample(1:length(sites), nrow(lung), replace=T))
#'######################### setup  ODAL #############################
#'control <- list(project_name = 'Lung cancer study',
#'        step = 'initialize',    # current step, updated by lead
#'        sites = sites,
#'        heterogeneity = FALSE,
#'        model = 'ODAL',
#'        outcome = "status",
#'        variables = c('age', 'sex'),
#'        optim_maxit=100,
#'        lead_site = sites[1],
#'        upload_date = as.character(Sys.time()),
#'        heterogeneity = FALSE)
#'## RUN BY LEAD ONLY 
#'pda(site_id='site1',control=control)
#'#run pda until step is empty
#'while (is.character(control$step)) {
#'  print(paste("step:",control$step))
#'  #cycle through sites
#'  for(i in 1:length(sites)) {
#'    control<-pda(ipdata=lung_split[[i]],site_id=sites[i])
#'  }
#'}
#'
#' @export
pda <- function(ipdata=NULL,site_id,control=NULL,dir=NULL,uri=NULL,secret=NULL){
  cloud_config<-getCloudConfig(site_id,dir,uri,secret)
  if(!(is.null(control)) &&  cloud_config$site_id==control$sites[1]) {
           pdaPut(control,'control',cloud_config)
           return(control)
  }
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
    #check if lead and sync needed
    if(cloud_config$site_id==control$sites[1]) {
       files<-pdaList(cloud_config) 
       if(all(paste0(control$sites,"_",control$step) %in% files)){
           control<-pdaSync(cloud_config)
       } 
    }
  }
  control
}


#' update the PDA control, used by the master site
#' @usage pdaSync <- function(config,site_id)
#' 
#' @param control_update Logical, update the PDA control?
#' 
pdaSync <- function(cloud_config){
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
  control
}
