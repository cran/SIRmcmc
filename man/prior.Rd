\name{prior}
\alias{prior}
\title{Compute prior probability of parameters}
\description{Compute the prior probability of alpha, beta, gamma, and epsilon}
\usage{
prior(alpha, beta, gamma, epsilon, esteps, params = list(alpha =
      list(location = 0, scale = 2), beta = list(shape = 0.01, rate = 0.01),
      gamma = list(shape = 0.01, rate = 0.01),
      epsilon = list(location = 0, scale = 2)))
}
\arguments{
  \item{alpha}{
    A number. The logarithm of alpha follows a logistic distribution.
  }
  \item{beta}{
    A number. beta follows a gammadistribution.
  }
  \item{gamma}{
    A number. gamma follows a gamma distribution.
  }
  \item{epsilon}{
    A number. The logarithm of epsilon follows a logistic distribution.
  }
  \item{esteps}{
    A number in (0,1). Used for logic as to whether to compute the
    prior probability of epsilon (1) or not (0).
  }
  \item{params}{
    A list of parameters for the prior distributions.
  }
}
\value{
  A probability of the model parameters under the prior distributions.
}
\seealso{
  \code{\link{dlogis}}
  \code{\link{dgamma}}
  \code{\link{household_transmission}}
}
