% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{LOS}
\alias{LOS}
\title{Length of Stay data}
\format{
A data frame with 1000 rows and 5 variables:
\describe{
  \item{site}{site id, 500 'site1', 400 'site2' and 100 'site3'}
  \item{age}{3 categories, 'young', 'middle', and 'old'}
  \item{sex}{2 categories, 'M' for male and 'F' for female}
  \item{lab}{lab test results, continuous value ranging from 0 to 100}
  \item{los}{LOS in days, ranging from 1 tp 28. Treated as continuous outcome in DLM}
}
}
\usage{
LOS
}
\description{
A simulated data set of hospitalization Length of Stay (LOS) from 3 sites
}
\keyword{datasets}
