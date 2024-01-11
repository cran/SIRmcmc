\name{household_transmission}
\alias{household_transmission}
\title{Estimate parameters from SIR model}
\description{Use the Metropolis algorithm to estimate parameters from
  the SIR compartmental model}
\usage{
household_transmission(onset_date, household_index, covariate = NULL,
                       followup_time, iterations,
                       delta = c(0.1,0.3,0.4,0.1), plot_chain = TRUE,
                       index = 1, start_date = NULL, prior = unifprior,
                       constant_hazard = FALSE,...)
}
\arguments{
  \item{onset_date}{
    A vector of onset dates. If onset_date is a character vector,
    then it will be coerced to dates using
    \code{\link{as.Date}}. \code{NA} values are interpretted as remaining
    susceptible at the end of followup.
  }
  \item{household_index}{
    A vector identifying households, which should be
    the same length as the onset_date argument.
  }
  \item{covariate}{
    A vector identifying covariate patterns. If given, then it is
    interpretted as a factor. A value for epsilon will be given for each
    level of covariate. If \code{NULL}, then epsilon is not estimated.
  }
  \item{followup_time}{
    An integer for the followup time in days.
  }
  \item{iterations}{
    An integer for the number of iterations of the Metropolis algorithm.
  }
  \item{delta}{
    A vector of length 4 for tuning the acceptance rate of the
    Metropolis algorithm. The order is (1) alpha, (2) beta, (3) gamma,
    and (4) epsilon.
  }
  \item{plot_chain}{
    A boolean of whether to plot the value of the chain vs the iterate.
  }
  \item{index}{
    An integer for the index date. Probabilities are conditional on the
    index date.  Any coordinates of onset_date equal
    to index will have a likelihood of 1. If you want unconditional
    probabilities, then index should be less than start_date.
  }
  \item{start_date}{
    Should be the same format as onset_date. Specifies the start of the
    followup period.
  }
  \item{prior}{
    A function to compute the prior probability of alpha, beta, gamma,
    and epsilon. Any user written function must take the arguments
    alpha, beta, gamma, epsilon, and esteps. Builtin functions are
    unifprior and \code{\link{prior}}. unifprior is uniform on (0.01,1000).
  }
  \item{constant_hazard}{
    If \code{FALSE}, then the algorithm computes a time dependent hazard
    for the cohort and the hazard from the community is proportional to
    the hazard in the cohort. If \code{TRUE}, then the algorithm assumes a
    constant hazard from the community.
  }
  \item{...}{
    Arguments to be passed to the function in the prior argument.
  }
}
\value{
  An object of the class \code{SIRmcmc}.
}
\details{
  If no covariates are supplied, only the model parameters alpha, beta,
  and gamma are estimated using a stepwise Metropolis algorithm. The
  model parameters are drawn from a uniform distribution on (0.1,1). The
  first step proposes a new alpha using \code{rnorm} and the mixing parameter
  \code{delta[1]}. The second and third steps are similar for beta and
  gamma. If covariates are supplied, the additional parameters,
  collectively called epsilon, are also estimated.

  The algorithm assumes followup begins on start_date and lasts for
  followup_time days. Any coordinate of onset_date equal to index does
  not contribute to the likelihood.

  Two priors are builtin: \code{\link{prior}} and
  \code{\link{unifprior}}. User defined prior functions must take the
  arguments alpha, beta, gamma, epslion, and esteps.
}
\examples{

##A trivial example-------------------------------------------------------
library(graphics)
onset<-sample(c(seq(1,10),rep(Inf,20)),size=500,replace=TRUE)
hh<-sample(seq(1,300),size=500,replace=TRUE)
chain<-household_transmission(onset_date = onset, household_index = hh,
                              followup_time = 10, iterations = 100)
community_attack_rate(SIRmcmc=chain)
secondary_attack_rate(household_size=3,SIRmcmc=chain)

##An example with household transmission---------------------------------
library(graphics)
data(hh.transmission)
set.seed(1)
iterations<-100
T<-30
delta<-c(0.1,0.6,0.8)
index<-0
##Find the MCMC estimates of alpha, beta, and gamma
chain<-household_transmission(
    onset_date=hh.transmission$onset,
    household_index=hh.transmission$household,
    covariate=NULL,
    followup_time=T,
    iterations=iterations,
    delta=delta,
    prior=unifprior,
    index=index
)
#Tabulate true type of transmission
hh_table<-table(
    table(
        is.finite(MCMC_date(hh.transmission$onset)),
        hh.transmission$household)["TRUE",]
)
##Calculate the true SAR
truth_table<-table(hh.transmission$transmission)
truth<-unname(truth_table["Household"]/sum(hh_table[2:3]))
cat("\n\nTrue Value of SAR\n\n")
print(truth)
##Find point and 95\% central creditable intervals for MCMC SAR
cat("\n\nMCMC Estimate of SAR\n\n")
secondary_attack_rate(household_size=2,SIRmcmc=chain)
days<-NULL
for(d in c(seq(1:5))){
    days<-c(days,as.character(d))
    a<-sum(table(tapply(X=hh.transmission$onset,INDEX=hh.transmission$household,FUN=diff))[days])
    cat(
        paste0(
            "\n\n",
            d,
            " Day Counting Estimate of SAR\n\n"
        )
    )
    ##Find point and 95\% confidence intervals for normal approx to SAR
    print(
        a/sum(hh_table[2:3])+c(p=0,LB=-1,UB=1) *
        qnorm(p=0.975) *
        sqrt(a*(hh_table[2]+hh_table[3]-a)/(hh_table[2]+hh_table[3])^3)
    )
}


##An example with rate ratios----------------------------------------
\dontrun{
    library(graphics)
    data(hh.transmission.epsilon)
    set.seed(1)
    iterations<-100
    T<-30
    delta<-c(0.1,0.1,0.1,0.1)
    index<-0
    ##Find the MCMC estimates of alpha, beta, gamma, and epsilon
    chain<-household_transmission(
        onset_date=hh.transmission.epsilon$onset,
        household_index=hh.transmission.epsilon$household,
        covariate=hh.transmission.epsilon$epsilon,
        followup_time=T,
        iterations=iterations,
        delta=delta,
        prior=unifprior,
        index=index
    )

    ##Find point and 95\% central creditable intervals for MCMC SAR
    cat("\n\nMCMC Estimate of SAR\n\n")
    print(secondary_attack_rate(household_size=2,SIRmcmc=chain))
    for(e in c(1,5)){
        hh_table<-table(table(
            is.finite(MCMC_date(hh.transmission.epsilon$onset)),
            hh.transmission.epsilon$household,
            hh.transmission.epsilon$epsilon)["TRUE",,as.character(e)])
        ##Tabulate true type of transmission
        truth_table<-table(
            hh.transmission.epsilon$transmission[which(
                hh.transmission.epsilon$epsilon==e
            )])
        ##Calculate the true SAR
        truth<-unname(truth_table["Household"]/sum(hh_table[2:3]))
        cat("\n\nTrue Value of SAR\n\n")
        print(truth)
        days<-NULL
        for(d in c(seq(1:5))){
            days<-c(days,as.character(d))
            a<-sum(table(tapply(
                X=hh.transmission.epsilon$onset[which(hh.transmission.epsilon$epsilon==e)],
                INDEX=hh.transmission.epsilon$household[which(hh.transmission.epsilon$epsilon==e)],
                FUN=diff))[days])
            ##Find point and 95\% confidence intervals for normal approx to SAR
            cat(paste0(
                "\n\n",
                d,
                " Day Counting Estimate of SAR\n\n"
            ))
            print(
                a/sum(hh_table[2:3])+c(p=0,LB=-1,UB=1)
                * qnorm(p=0.975)
                * sqrt(a*(hh_table[2]+hh_table[3]-a)/(hh_table[2]+hh_table[3])^3)
            )
        }
    }
}

}

\seealso{
  \code{\link{prior}}
  \code{\link{unifprior}}
}