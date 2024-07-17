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

# set in pda() ?
dGEM.steps <- c('initialize','derive','estimate','synthesize')

# Note: The implementation of downstream analysis (i..e, hospital profiling) will 
# depend on the guidelines specified in the protocol, 
# allowing for the optional implementation in alignment with project-specific needs. 
# If a project's emphasis is on deploying a decentralized GLMM without hospital profiling, 
# there is no necessity to proceed to 'estimate' step, as after 'derive', the coordianting
# center will have the esimated common effects and random effects. 
dGEM.family <- 'binomial'

#' @useDynLib pda
#' @title dGEM initialize
#' 
#' @usage dGEM.initialize(ipdata,control,config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#' 
#' @references NA
#' @return init
#' @keywords internal
dGEM.initialize <- function(ipdata,control,config){
  fit_i <- glm(status ~ 0+., data=ipdata,family = "binomial"(link = "logit"))  
  fit_null <- glm(status ~ 1, data = ipdata, family = 'binomial')
  R2 <- 1 - logLik(fit_i)/logLik(fit_null)
  AIC <- AIC(fit_i)
  BIC <- BIC(fit_i)


  ### Split the data and then estimate the AUC
  # Number of repetitions
  n_repeats <- 50
  
  # Initialize a vector to store AUC values
  auc_values <- numeric(n_repeats)
  
  # Repeat the process 50 times
  set.seed(1234) # Set seed for reproducibility
  for (i in 1:n_repeats) {
    # Split the data into training and testing
    sample_size <- floor(0.75 * nrow(ipdata))
    train_indices <- sample(seq_len(nrow(ipdata)), size = sample_size)
    train_data <- tmp_data_to_fit[train_indices, ]
    test_data <- tmp_data_to_fit[-train_indices, ]
    
    # Check if test_data$outcome has two levels
    if (length(unique(test_data$status)) < 2) {
      next # Skip this iteration
    }
    
    # Fit logistic regression model on training data
    model <- glm(status ~ ., data = train_data, family = binomial)
    
    # Make predictions on the testing data
    predictions <- predict(model, newdata = test_data, type = "response")
    
    # Calculate AUC
    roc_curve <- roc(test_data$status, predictions)
    auc_value <- auc(roc_curve)
    
    # Store AUC value
    auc_values[i] <- auc_value
  }
  
  # Calculate the average AUC
  average_auc <- mean(auc_values)
  
  init <- list(site = config$site_id,
               site_size = nrow(ipdata),
               bhat_i = fit_i$coef,
               Vhat_i = diag(vcov(fit_i)),
              R2 = R2,
              auc = average_auc,
              AIC = AIC,
              BIC = BIC)  # glm summary(fit_i)$coef[,2]^2 may omit NA's
  return(init)
}

#' @useDynLib pda
#' @title dGEM hospital-specific effect derivation
#' 
#' @usage dGEM.derive(ipdata,control,config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#' @param hosdata hospital-level data
#'
#' @return  hospital_effect
#' @keywords internal
dGEM.derive <- function(ipdata,control,config,hosdata){
  
  # b_meta <- betameta
  bbar <- control$beta_init[-1]
  
  # calculate the site-specific effects with beta as offset
  os = as.matrix(ipdata[,-c(1,2)]) %*% as.matrix(bbar)
  fit_i = glm(status ~ 1 + offset(os),
              data = ipdata[,-2], family = "binomial")
  fit_summary = summary(fit_i)
  
  hospital_effect <- list(
    site=config$site_id, 
    site_size = nrow(ipdata),
    hosdata = hosdata,
    gammahat_i=fit_summary$coefficients[1],
    Vgammahat_i=fit_summary$coefficients[2]^2)
  
  return(hospital_effect)
}


#' @useDynLib pda
#' @title dGEM standardized event rate estimation
#' 
#' @usage dGEM.estimate(ipdata,control,config)
#' @param ipdata local data in data frame
#' @param control PDA control
#' @param config cloud configuration
#' 
#' @details step-3: 
#' @return  event rate
#' @keywords internal
dGEM.estimate <- function(ipdata,control,config) {
  
  hospital_effect <- pdaGet('control',config)$estimated_hospital_effect
  bbar <- control$beta_init[-1]
  
  # calculate event rate
  event_rate = c()
  ind = 1
  for (site_i in control$sites){
    gamma_i <- hospital_effect[ind]
    val = gamma_i + as.matrix(ipdata[,-c(1,2)]) %*% as.matrix(bbar)
    tmp = exp(val)/(1+exp(val))
    event_rate[ind] = sum(tmp)
    ind = ind + 1
  }
  
  
  event_rate <- list(sum_event_rate = event_rate,
                     site=config$site_id, 
                     site_size=nrow(ipdata))
  ######################################################
  
  return(event_rate)
}



#' @useDynLib pda
#' @title PDA dGEM synthesize 
#' 
#' @usage dGEM.synthesize(ipdata,control,config)
#' @param control pda control
#' @param config pda cloud configuration
#' 
#' @details Synthesis to get the standardized mortality rate 
#'
#' @return  list(final_event_rate=final_event_rate)
#' @keywords internal
dGEM.synthesize <- function(control,config) {
  
  N = 0
  for (site_i in control$sites){
    n =  pdaGet(paste0(site_i,'_estimate'),config)$site_size
    N = N + n
  }
  
  event_rate_mat = pdaGet(paste0(control$sites[1],'_estimate'),config)$sum_event_rate
  for (site_i in control$sites[-1]){
    event_rate_mat =  rbind(event_rate_mat, pdaGet(paste0(site_i,'_estimate'),config)$sum_event_rate)
  }
  
  
  # Final standardized event rate
  final_event_rate = apply(event_rate_mat, 2, sum)/N
  names(final_event_rate) = control$sites
  
  return(list(final_event_rate=as.list(final_event_rate)))
}
