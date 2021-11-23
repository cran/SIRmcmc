\name{hh.transmission.epsilon}
\alias{hh.transmission.epsilon}
\docType{data}
\title{
  Simulated data under SIR, with covariates.
}
\description{
 Simulated data under the SIR model with a high risk and low risk
 group. The first 1000 observations are identical to the data
 \code{hh.transmission}.
}
\usage{data("hh.transmission.epsilon")}
\format{
  A data frame with 2000 observations on the following 4 variables.
  \describe{
    \item{\code{household}}{A numeric vector of household IDs.}
    \item{\code{onset}}{Date of onset for cases.}
    \item{\code{transmission}}{A factor with levels \code{Community}
      \code{Household} \code{None} representing source of infection.}
    \item{\code{epsilon}}{A numeric vector of the value of epsilon used
      to simulate the data.}
  }
}