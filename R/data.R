
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