\name{MCMC.date}
\alias{MCMC.date}
\title{Convert dates to a list of extendend natural numbers.}
\description{Converts dates to a list of numbers representing the number
of days from the start of followup until the start of the infectious period.}
\usage{
MCMC.date(dates,start.date=NULL)
}
\arguments{
  \item{dates}{
    A list of dates. May be of character, numeric, or date class.
  }
  \item{start.date}{
    The start date of follow up
  }
}
\details{
  Covrerts dates to days of followup until the start of the infectious
  period. Missing data are set to infinity and are assumed susceptible
  during followup.
}
\value{
  A list of extended natural numbers.
}