\name{hh.transmission}
\alias{hh.transmission}
\docType{data}
\title{
Simulated data under the SIR.
}
\description{
Simulated data under the SIR model, without covariates.
}
\usage{data("hh.transmission")}
\format{
  A data frame with 1000 observations on the following 3 variables.
  \describe{
    \item{\code{household}}{A numeric vector of household IDs.}
    \item{\code{onset}}{Date of onset for cases}
    \item{\code{transmission}}{A factor with levels \code{Community}
      \code{Household} \code{None} representing source of infection.}
  }
}