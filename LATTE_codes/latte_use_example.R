source("/Users/luli/Downloads/LATTE_final_results/001_lossless_property/0_rscripts/helper_iptw.r")
source("/Users/luli/Downloads/LATTE_final_results/001_lossless_property/0_rscripts/helper_strat.R")





drugs_path <- "/Volumes/home_luli1/fed_em/cleaned_codes_for_other_sites/dementia_upenn_min/outcomes_all_100_now_upenn/save_cohort_all_loose2_5yr_with_index_race/cohort_all_name_size_positive.csv"
drugs <- read.csv(drugs_path)
drugs <- drugs[drugs$n_patients >= 100, ]
ids <- drugs[, 2]
 
drug_names = c()
latte_ests = c()
latte_ses = c()
normal_ests = c()
normal_ses = c()

for (id in ids){
drug_name = drugs[drugs$cohort_name == id, "drug_name"]
if (!drug_name %in% c("donepezil")) {
    id = str_replace(id, ".pkl", "")
    drug_names = c(drug_names, drug_name)
    data_path = paste0("/Volumes/home_luli1/fed_em/cleaned_codes_for_other_sites/dementia_upenn_min/results_alldata_5yr_w_index_race_now_upenn_supplement_100/", id, "all_data.csv")
    cohort = read.csv(data_path)
    cohort = merge_duplicates(cohort)
    cohort <- cohort %>% select(-X)


    cohort <- cohort %>%
        mutate(across(ends_with("_value"), ~ replace(., . == -1, 0)))

    ##########################################################################
    #########               Propensity Score Model            ################
    ##########################################################################
    xvars <- colnames(cohort)[!grepl("^outcome_", colnames(cohort)) & !colnames(cohort) %in% c("ID", "treatment", "index_date")]
    xvars = xvars[colSums(cohort[xvars]) > 30]
    yvars <- colnames(cohort)[grepl("^outcome_", colnames(cohort))]
    form <- as.formula(paste("treatment ~ ", paste(xvars, collapse = "+"))) # model formula
    mydata_ps <- cohort[, colnames(cohort) %in% c(xvars, "treatment", yvars)]
    Xmat <- grab_design_matrix(data = mydata_ps, rhs_formula = form)
    Y <- mydata_ps$treatment

    # Fit a GLM
    nfolds <- 10
    foldid <- sample(rep(seq(nfolds), length.out = length(Y)))
    Fit_ps_cv <- cv.glmnet(Xmat, Y, alpha = 1, family = "binomial", nfolds = nfolds, foldid = foldid) # cross validation to select lambda
    Fit_ps <- glmnet(Xmat, Y, alpha = 1, family = "binomial", lambda = Fit_ps_cv$lambda.min)
    propensityScore <- predict(Fit_ps, (Xmat), type = "response")
    mydata_ps$propensityScore <- propensityScore


    best_n_strata = 0
    best_after = 1000
    best_smd_res = NULL
    for (nstrata in 2:6) {
        stratifiedPop = get_stratified_pop(mydata_test = mydata_ps, nstrata = nstrata)
        smd_res <- get_SMD(stratifiedPop = stratifiedPop, xvars = xvars)
        before = sum(abs(smd_res$smd_before$SMD) > 0.2)
        after = sum(abs(smd_res$smd_after$SMD) > 0.2)
        if (after < best_after) {
            best_after = after
            best_n_strata = nstrata
            best_smd_res = smd_res
        }
    }


    stratifiedPop = get_stratified_pop(mydata_test = mydata_ps, nstrata = best_n_strata)




    ###### LATTE results
    outcome_id = "ADRD"

    outcome_formula <- as.formula(paste0("outcome_", outcome_id, "_value ~ treatment"))
    ADdata = NULL
    ADdata <- create_2x2_tables(stratifiedPop$treatment, stratifiedPop[[paste0("outcome_", outcome_id, "_value")]], stratifiedPop$stratumId)


    outcome_formula <- as.formula("Y ~ treatment")
    KSiteAD_uf <- ADdata
    ResAD <- optimize_conditional_logistic_2x2(KSiteAD_uf)
    confint <- ResAD$coefficients + c(-1, 1) * 1.96 * ResAD$se
    # Store results in a list
    results <- list(
        drug_id = paste0("435", " (", "435", ")"),
        se = (ResAD$se),
        est = (ResAD$coefficients),
        ll = (confint[1]),
        ul = (confint[2])
    )
    latte_est = results$est
    latte_se = results$se

    outcome_formula = "outcome_ADRD_value ~ treatment + stratumId"
    ##### directly fitting stratified losgitic regression
    fit <- glm(outcome_formula, data = stratifiedPop, family = binomial(link = "logit"))
    normal_est = summary(fit)$coefficients[2, 1]
    normal_se = summary(fit)$coefficients[2, 2]

    # plot the results
    latte_ests = c(latte_ests, latte_est)
    latte_ses = c(latte_ses, latte_se)
    normal_ests = c(normal_ests, normal_est)
    normal_ses = c(normal_ses, normal_se)
}
}

### save the results
results_df <- data.frame(
    drug_id = drug_names,
    latte_ests = latte_ests,
    latte_ses = latte_ses,
    normal_ests = normal_ests,
    normal_ses = normal_ses
)

write.csv(results_df, "method_comparison.csv")
original_drug_names = drug_names
original_latte_ests = latte_ests
original_latte_ses = latte_ses
original_normal_ests = normal_ests
original_normal_ses = normal_ses

indices = which(original_drug_names %in% c("gabapentin", "omeprazole", "pantoprazole", "famotidine", "methylprednisolone", "amlodipine"))
drug_names = original_drug_names[indices]
latte_ests = original_latte_ests[indices]
latte_ses = original_latte_ses[indices]
normal_ests = original_normal_ests[indices]
normal_ses = original_normal_ses[indices]

# Helper function to optimize number of strata
optimize_strata <- function(data, xvars, min_strata = 2, max_strata = 6) {
  best_n_strata <- 0
  best_after <- 1000
  best_smd_res <- NULL
  
  for (nstrata in min_strata:max_strata) {
    stratifiedPop <- get_stratified_pop(data, nstrata = nstrata)
    smd_res <- get_SMD(stratifiedPop = stratifiedPop, xvars = xvars)
    
    before <- sum(abs(smd_res$smd_before$SMD) > 0.2)
    after <- sum(abs(smd_res$smd_after$SMD) > 0.2)
    
    if (after < best_after) {
      best_after <- after
      best_n_strata <- nstrata
      best_smd_res <- smd_res
    }
  }
  
  return(list(
    n_strata = best_n_strata,
    smd_res = best_smd_res
  ))
}
