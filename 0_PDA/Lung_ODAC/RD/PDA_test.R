


require(survival)
require(data.table)
require(Rcpp)
# require(ODACO)

# rm(list=ls())


################################ this session requires your input ###################################################
mysite <- 'CL'    # provide your name initial
mycloud <- '/Users/chongliangluo/Dropbox/0_PDA/Lung_ODAC/cloud/'     # provide your Dropbox PDA cloud folder 
#####################################################################################################################

setwd(paste0(mycloud, '../', mysite))

Rcpp::sourceCpp('../../PDA/src/rcpp_coxph.cpp')
source('../../PDA/R/PDA_engine.R')


## prepare test data at your site
data(lung)
lung2=lung[,1:5]
lung2$inst[lung$inst>=2] = 2
lung2$inst[lung$inst>=12] = 3
lung2$inst[is.na(lung2$inst)] <- 3
lung2$inst <- c('CL', 'RD', 'JT')[lung2$inst]   # LETTERS[lung2$inst]
table(lung2$inst)
lung2$status <- lung2$status == 2
lung2$sex <- lung2$sex-1
lung2 <- data.table(lung2)
fit.pool <- coxph(Surv(time, status) ~ age+sex, data=lung2)
fit.pool$coef

mydata <- lung2[inst==mysite, -'inst']


## research project setting, decided by all collaborators 
pda_control <- list(cloud = mycloud,       # specify the cloud folder at your PC
                    outcome = 'Lung',
                    risk_factor = c('age', 'sex'),
                    model = 'ODAC', 
                    all_site = c('CL', 'RD', 'JT'),
                    local_site = mysite,   # specify your site as local
                    optim_maxit=100)

# pda_control$cloud <- mycloud           # specify your cloud folder  
# pda_control$local_site <- mysite       # specify your site as local



# step from
# c(1:5) or
# c('initialize', 'summary_stat', 'derivatives', 'surrogate_est', 'synthesize')

pda_main(mydata, step = 1, control = pda_control)

# wait until all sites completed step-1
pda_main(mydata, step = 2, control = pda_control)


# wait until all sites completed step-2
pda_main(mydata, step = 3, control = pda_control)


# wait until all sites completed step-3
output <- pda_main(mydata, step = 4, control = pda_control)
pda_broadcast(output, obj_type = 'surrogate_est', control = pda_control)


# wait until all sites completed step-4 (optional)
output <- pda_main(mydata, step = 5, control = pda_control)





















# fit.ODAC <- DistCox(mydata = lung2, id.local = mysite, init_est = b_meta, strat = F)
# fit.local <- coxph(Surv(time, status) ~ age+sex, data=mydata)
# fit.pool <- coxph(Surv(time, status) ~ age+sex, data=lung2)








# ## 
# step <- 4
# for(site in pda_control$all_site){
# 
#   mydata <- lung2[inst==site, -'inst']
#   
#   pda_control <- list(cloud = mycloud,       # specify the cloud folder at your PC
#                       outcome = 'Lung',
#                       risk_factor = c('age', 'sex'),
#                       model = 'ODAC', 
#                       all_site = c('CL', 'RD', 'JT'),
#                       local_site = mysite,   # specify your site as local
#                       optim_maxit=100)
#   pda_control$local_site <- site
# 
#   if(step==1){
#     tmp <- pda_initialize(mydata, broadcast=TRUE, control=pda_control)
#     print(tmp$bhat_i)
#     # age         sex
#     # 0.04028874 -0.57977372
#     # age        sex
#     # 0.0157372 -0.1581604
#     # age          sex
#     # 0.006101769 -0.796974901
#   }
# 
#   if(step==2){
#     tmp <- pda_summary_stat(bbar=NULL, mydata, broadcast=TRUE, control=pda_control)
#     print(tmp$b_meta)
#     # 0.02029786 -0.63152898
#   }
# 
#   if(step==3) tmp <- pda_derivatives(bbar=NULL, mydata, broadcast=TRUE, control=pda_control)
# 
# 
#   if(step==4){
#     tmp <- pda_surrogate_est(mydata, bbar=NULL, broadcast=TRUE, control=pda_control)
#     print(tmp$btilde)
#     # age         sex
#     # 0.01697928 -0.51267777
#     # age         sex
#     # 0.01696673 -0.51262126
#     # age         sex
#     # 0.01703605 -0.51276307
#   }
# 
#   if(step==5) tmp <- pda_synthesize(control=pda_control)
# 
# }
# # }
#  

# compare: DistCox (beta_tilde and beta_tilde1) obtains coef est better than using local only, closer to use pooled data
# rbind(bhat_i=fit.local$coef,
#       btilde_i=fit.ODAC$beta_tilde,
#       # btilde_1i=fit.ODAC$beta_tilde1,
#       bpool=fit.pool$coef,
#       b_meta,
#       btilde_1=c(0.02640378, -0.58091989)
#       , btilde=c(tmp$btilde)
#       )

# bhat_i   0.04028874 -0.5797737
# btilde_i 0.01703068 -0.5131145
# bpool    0.01704533 -0.5132185
# b_meta   0.02029786 -0.6315290
# btilde_1 0.01697928 -0.51267777    # 0.02640378 -0.5809199
# btilde   0.01699429 -0.51269058    # 0.02636965 -0.5819797
 
