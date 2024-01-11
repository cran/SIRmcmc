\name{community_attack_rate}
\alias{community_attack_rate}
\title{Compute the community attack rate}
\description{Computes the community attack rate for a a cohort using the
  value of a call to \code{\link{household_transmission}}.
}
\usage{
community_attack_rate(SIRmcmc,probs=c(0.5,0.025,0.975))
}
\arguments{
  \item{SIRmcmc}{
    The value of a call to \code{\link{household_transmission}}.
  }
  \item{probs}{
    A numeric vector of the quantiles of the posterior distribution. The
    default is the median and the central 95\% credible region.
  }
}
\details{
  Computes the posterior probability distribution of the community
  attack rate from the Metropolis algoritm. Returns quantiles of the
  distribution specified in the \code{probs} argument.
}
\value{
  An array of community attack rates. The first dimension is the
  value of epsilon in \code{SIRmcmc}. The second dimension
  is the posterior probability in the \code{probs} argument
}