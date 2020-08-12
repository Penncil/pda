require(data.table)
require(Rcpp)
## master site update pda_control$step = 3
pda_control_update()
## surrogate_est
b_surr = pda(data = mydata, mysite='site1')
## compare
lung2 = fread('data/Lung.csv')
fit.pool <- glm(status ~ age + sex, data = lung2, family = "binomial")
fit.pool$coef
rbind(pooled= fit.pool$coef,
      meta = pda_control$beta_init,
      odal = b_surr$btilde) 
