\name{unifprior}
\alias{unifprior}
\title{A uniform prior on the model parameters}
\description{Compute a uniform prior on alpha, beta, gamma, and epsilon}
\usage{
unifprior(alpha, beta, gamma, epsilon, esteps, UB = 1000, LB = 0.01)
}
\arguments{
  \item{alpha}{
    A number compared to the upper bound UB and the lower bound LB.
  }
  \item{beta}{
    A number compared to the upper bound UB and the lower bound LB.
  }
  \item{gamma}{
    A number compared to the upper bound UB and the lower bound LB.
  }
  \item{epsilon}{
    A number compared to the upper bound UB and the lower bound LB.
  }
  \item{esteps}{
    A number in (0,1). Used for logic as to whether to compute the
    prior probability of epsilon (1) or not (0).
  }
  \item{UB}{
    A number used as an upper bound on the model parameters.
  }
  \item{LB}{
    A number used as a lower bound on the model parameters.
  }
}
\value{
  If all the model parameters are between the lower bound and the upper
  bound, then unifprior returns 1. Otherwise, unifprior returns 0.
}
\seealso{
  \code{\link{dlogis}}
  \code{\link{dgamma}}
  \code{\link{household_transmission}}
  \code{\link{prior}}
}
