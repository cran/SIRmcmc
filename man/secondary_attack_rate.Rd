\name{secondary_attack_rate}
\alias{secondary_attack_rate}
\title{Compute the secondary attack rate}
\description{Using the value of a call to
  \code{\link{household_transmission}}, computes the secondary attack
  rate for households.}
\usage{
secondary_attack_rate(household_size,SIRmcmc,probs=c(0.5,0.025,0.975))
}
\arguments{
  \item{household_size}{
    A numeric vector containing the number of people in a household.
  }
  \item{SIRmcmc}{
    The value of a call to \code{\link{household_transmission}}.
  }
  \item{probs}{
    A numeric vector of the quantiles of the posterior distribution. The
    default is the median and the central 95\% credible region.
  }
}
\details{
  Computes the posterior probability distribution of the secondary
  attack rate from the Metropolis algorithm. Returns quantiles of the
  distribution specified in the \code{probs} argument.
}
\value{
  An array with secondary attack rates.

  The first dimension of the array is the household size.

  The second dimension of the array is the quantiles of the posterior
  distribution in the \code{probs} argument.

  The third dimension of the array is the value of epsilon in \code{SIRmcmc}.
}