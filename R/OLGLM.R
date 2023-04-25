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
OLGLM.steps <- c('initialize')
OLGLM.family <- 'binomial'

#' @useDynLib pda
#' @title ODAL initialize
#' 
#' @usage ODAL.initialize(ipdata,control,config)
#' @param ipdata individual participant data
#' @param control pda control data
#' @param config local site configuration
#' 
#' @references 
#' @return init
#' @keywords internal
OLGLM.initialize <- function(ipdata,control,config){

  
  Xmat <- ipdata[,!colnames(ipdata) %in% control$outcome, with = FALSE] 
  Y <- ipdata[,colnames(ipdata) %in% control$outcome, with = FALSE]
  Xmat.tbl <- data.frame(Xmat)
  category_combinations <- expand.grid(lapply(Xmat.tbl, unique),
                                       stringsAsFactors = FALSE)%>%arrange_all()
  colnames(Xmat.tbl) <- colnames(category_combinations)
  
  if(control$link == "canonical"){
    Xmat.tbl <- as_tibble(Xmat.tbl)
    cols <- colnames(Xmat.tbl)
    Xtable_initial <- Xmat.tbl %>%group_by_at(.vars = cols)%>%summarise(n = n())
    Xtable <- category_combinations%>%left_join(Xtable_initial,by = cols) %>% as.data.frame()
    Xtable$n[which(is.na(Xtable$n))] = 0
    colnames(Xtable) <- c(colnames(Xmat),'n')
    if(is.numeric(control$cutoff)==TRUE){
      Xtable$n[Xtable$n>0 & Xtable$n <control$cutoff] <- rep(ceiling(control$cutoff/2))
    }
    SXY <- t(Y)%*%as.matrix(Xmat) 
    AD <- list(SXY = SXY, Xtable = Xtable)
  }
  return(AD)
}

