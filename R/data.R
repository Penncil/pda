
#' Lung cancer survival time data
#'
#' A data set modified from the lung data in survival package (see demo(ODAC)).
#'
#' @format A data frame with 228 rows and 5 variables:
#' \describe{
#'   \item{site}{simulated site id, 86 'site1', 83 'site2' and 59 'site3'}
#'   \item{time}{survival time in days}
#'   \item{status}{censoring status 0=censored, 1=dead}
#'   \item{age}{age in years}
#'   \item{sex}{1 for female and 0 for male}
#' }
#' @source \url{https://CRAN.R-project.org/package=survival}
"lung2"


#' CrabSatellites data
#'
#' A data set modified from the CrabSatellites data in countreg package (see demo(ODAH)).
#'
#' @format A data frame containing 173 observations on 4 variables.
#' \describe{
#'   \item{site}{Simulated site id, 85 'site1' and 88 'site2'.}
#'   \item{satellites}{Number of satellites. Treated as (zero-inflated) count outcome in ODAH}
#'   \item{width}{Carapace width (cm).}
#'   \item{weight}{Weight (kg).} 
#' }
#' @source \url{https://rdrr.io/rforge/countreg/man/CrabSatellites.html}
"cs"


#' Length of Stay data
#'
#' A simulated data set of hospitalization Length of Stay (LOS) from 3 sites
#'
#' @format A data frame with 1000 rows and 5 variables:
#' \describe{
#'   \item{site}{site id, 500 'site1', 400 'site2' and 100 'site3'}
#'   \item{age}{3 categories, 'young', 'middle', and 'old'}
#'   \item{sex}{2 categories, 'M' for male and 'F' for female}
#'   \item{lab}{lab test results, continuous value ranging from 0 to 100}
#'   \item{los}{LOS in days, ranging from 1 tp 28. Treated as continuous outcome in DLM}
#' } 
"LOS"


#' COVID-19 LOS and mortality data
#'
#' A simulated data set of hospitalization Length of Stay (LOS) and mortality from 6 sites
#'
#' @format A data frame with 2100 rows and 6 variables:
#' \describe{
#'   \item{site}{site id, 600 'site1', 500 'site2', 400 'site3', 300 'site4', 200 'site5', 100 'site6'}
#'   \item{age}{continuous age in year, min 3 max 97}
#'   \item{sex}{2 categories, '1' for male and '0' for female}
#'   \item{lab}{lab test results, continuous value ranging from 2.3 to 97.4}
#'   \item{los}{LOS in days, ranging from 1 to 29}
#'   \item{death}{mortality status, '1' for death and '0' for alive.}
#' } 
"covid"


#' ADAP simulated data 
#'
#' A simulated data set for ADAP demonstration
#'
#' @format A list containing the following elements:
#' \describe{
#'   \item{sites}{site id, 300 'site1', 300 'site2', 300 'site3'}
#'   \item{status}{binary outcome of length 900}
#'   \item{x}{900 by 49 matrix generated by standard normal distribution, representing the covariates} 
#' } 
"ADAP_data"

#' ODACAT simulated data 
#'
#' A simulated data set for ODACAT demonstration
#'
#' @format A data frame with 300 rows and 5 variables:
#' \describe{
#'   \item{id.site}{site id, 105 'site1', 105 'site2', 90 'site3'}
#'   \item{outcome}{3-category outcome, possible values are 1,2,3. Category 3 will be used as reference}
#'   \item{X1}{the first covariate, continuous}
#'   \item{X2}{the second covariate, binary}
#'   \item{X3}{the third covariate, binary}
#' } 
"ODACAT_ordinal"

#' ODACH_CC simulated data 
#'
#' A simulated data set for ODACH with case-cohort design demonstration
#'
#' @format A data frame with 413 rows and 8 variables:
#' \describe{
#'   \item{site}{site id, 187 'site1', 133 'site2', 93 'site3'. The full_cohort_size are 800, 600 and 400 respectively} 
#'   \item{subcohort}{1=subcohort, e.g. uniformly subsampled from full_cohort_size, 0=case} 
#'   \item{time}{survival time}
#'   \item{status}{censoring status 0=censored, 1=dead} 
#'   \item{X1}{the first covariate, continuous}
#'   \item{X2}{the second covariate, continuous} 
#'   \item{Category}{the third covariate, categorical} 
#'   \item{Group}{the fourth covariate, categorical} 
#' } 
"odach_cc"