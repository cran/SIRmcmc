##Try to coerce dates to an integer format for the MCMC algorithm
MCMC_date<-function(dates,start_date=NULL){
    if(is(dates, "character")){
        if(all(is.na(dates)==is.na(as.Date(dates)))==FALSE){
            failed_dates<-dates[is.na(dates)!=is.na(as.Date(dates))]
            failed_dates<-unique(failed_dates)
            stop(paste(
                "Coercing dates from character format to date format failed for the following:\n\n",
                failed_dates,
                "\n\nTry using as.Date() on your data."
            ))
        }else{
            dates<-as.Date(dates)
        }
    }
    if(is.null(start_date)==TRUE){
        dates<-as.numeric(dates-min(dates,na.rm=TRUE))+1
    }else{
        dates<-as.numeric(dates-start_date+1)
    }
    if(any(is.na(dates))==TRUE){
        warning("Assuming NA are non-cases")
        dates[which(is.na(dates))]<-Inf
    }
    return(dates)
}

##A lookup table for endpoints for integrating. Avoids redundant calls to unique().
make_cut_points<-function(onset_date,household_index,followup_time){
    cut_points<-list()
    length(cut_points)<-max(household_index)
    for(h in household_index){
        index<-which(household_index==h)
        temp_onset_date<-onset_date[index]
        if(all(is.finite(temp_onset_date))==FALSE){
            cut_points[[h]]<-sort(
                unique(
                    c(
                        0,
                        temp_onset_date,
                        temp_onset_date-1,
                        followup_time
                    )
                )
            )
        }else{
            cut_points[[h]]<-sort(
                unique(
                    c(
                        0,
                        temp_onset_date,
                        temp_onset_date-1
                    )
                )
            )
        }
    }
    return(cut_points)
}


##Computes cummulative hazard function within cohort. Returns a smoothed cummulative hazard and
##the product limit estimates of the cummulative hazard
make_cohort_hazard_table<-function(onset_date,followup_time,constant_hazard=FALSE){
    if(any(onset_date[which(is.finite(onset_date))]>followup_time)==TRUE){
        stop(paste("Onset date after followup ends. Rerun with followup >=",max(onset_date[which(is.finite(onset_date))])))
    }
    if(all(is.finite(onset_date))==FALSE){
        cut_points<-c(seq(0,followup_time),Inf)
    }else{
        if(all(onset_date<followup_time)==TRUE){
            stop(paste("100% attack rate before end of followup. Tried to evaluate Inf-Inf. Rerun with followup =",max(onset_date)))
        }
        cut_points<-c(seq(0,followup_time))
    }
    K<-length(cut_points)
    cases<-table(
        factor(
            x=onset_date,
            levels=cut_points
        )
    )
    survival<-rep(1,K)
    for(k in 2:K){
        survival[k]<-survival[k-1]*(1-cases[k]/sum(cases[k:K]))
    }
    Lambda<--log(survival)
    lambda<-diff(Lambda)
    if(length(K)<=7){
        if(is(try(loess(lambda[1:(K-2)]~cut_points[1:(K-2)])), "try-error")){
            smooth_lambda<-loess(lambda[1:(K-2)]~cut_points[1:(K-2)])$fitted
            if(any(smooth_lambda<0)){
                smooth_lambda<-loess(lambda[1:(K-2)]~cut_points[1:(K-2)],degree=0)$fitted
                warning("Using degree 0 loess")
            }
            if(any(smooth_lambda<0)){
                smooth_lambda<-lambda
                warning("Using empiric survival")
            }
        }else{
            smooth_lambda<-lambda
            warning("LOESS failed. Using empiric suvival")
        }
    }else{
        smooth_lambda<-lambda
        warning("Too few time points. Using empiric survival.")
    }
    smooth_Lambda<-cumsum(smooth_lambda)
    if(constant_hazard==FALSE){
        return(
            cbind(
                t=cut_points,
                Lambda=c(0,smooth_Lambda),
                empiric=Lambda
            )
        )
    }else{
        return(
            cbind(
                t=cut_points,
                Lambda=cut_points,
                empiric=Lambda
            )
        )
    }
}

##Sanity checks on log probabilities
check_probs<-function(p,q){
    if(is.na(p)==TRUE | is.na(q)){stop("Probability NA")}
    if(is.finite(p)==FALSE & is.finite(q)==FALSE){stop("Tried to Evaluate -Inf+Inf")}
    if(is.numeric(p)==FALSE | is.numeric(q)==FALSE){stop("Probability not a numeric variable")}
    if(is.nan(p)==TRUE | is.nan(q)==TRUE){stop(paste("Probability not a number"))}
}

##Find the index that sorts households, smallest households first.
sort_order_households<-function(household_index){
    n<-length(household_index)
    out_index<-rep(NA,length(n))
    hh_table<-sort(table(household_index))
    i<-1
    for(h in as.numeric(names(hh_table))){
        temp_index<-which(household_index==h)
        m<-length(temp_index)
        out_index[i:(i+m-1)]<-temp_index
        i<-i+m
    }
    return(out_index)
}

##Renumber the households in increasing order
renumber_households<-function(household_index){
    n<-length(household_index)
    out_index<-rep(NA,length(n))
    new_index<-1
    out_index[1]<-new_index
    for(i in 2:n){
        if(household_index[i]!=household_index[i-1]){
            new_index<-new_index+1
        }
        out_index[i]<-new_index
    }
    return(out_index)
}

##By default, a vague prior.
prior<-function(
                alpha,
                beta,
                gamma,
                epsilon,
                esteps,
                params=list(
                    alpha=list(location=0,scale=2),
                    beta=list(shape=0.01,rate=0.01),
                    gamma=list(shape=0.01,rate=0.01),
                    epsilon=list(location=0,scale=2)
                )){
    epsilon<-unique(epsilon)
    prob<-prod(
        dlogis(x=log(alpha),location=params$alpha$location,scale=params$alpha$scale),
        dgamma(x=beta,shape=params$beta$shape,rate=params$beta$rate),
        dgamma(x=gamma,shape=params$gamma$shape,rate=params$gamma$rate)
    )
    if(esteps==1){
        prob<-prod(
            prob,
            dlogis(x=log(epsilon),location=params$epsilon$location,scale=params$epsilon$scale)
        )
    }
    if(is.numeric(prob)==FALSE || is.nan(prob)){
        stop(
            paste(
                "Error in Prior: alpha =",
                alpha,
                "beta =",
                beta,
                "gamma =",
                gamma,
                "epsilon =",
                epsilon,
                "\n\nEvaluated to",
                prob
            )
        )
    }
    return(prob)
}

                                        #Metropolis algorithm
household_transmission<-function(
                                 onset_date,
                                 household_index,
                                 covariate=NULL,
                                 followup_time,
                                 iterations,
                                 delta=c(0.1,0.3,0.4,0.1),
                                 plot_chain=TRUE,
                                 index=1,
                                 start_date=NULL,
                                 prior=unifprior,
                                 constant_hazard=FALSE,
                                 ...
                                 ){
    ##Preprares the data for the algorithm: households together, small households first, and an increasing index.
    household_index<-as.integer(as.factor(household_index))
    onset_date<-MCMC_date(onset_date,start_date)
    sort_order<-sort_order_households(household_index)
    onset_date<-onset_date[sort_order]
    household_index<-household_index[sort_order]
    household_index<-renumber_households(household_index)
    ##Decide whether to estimate epsilon, and structure the covariate data for the algorithm.
    if(is.null(covariate)){
        covariate<-rep(1,length(onset_date))
        estimate_epsilon<-0
    }else{
        estimate_epsilon<-1
        covariate<-covariate[sort_order]
    }
    levels_covariate<-sort(unique(covariate))
    if(length(levels_covariate)==1 & estimate_epsilon==1){
        warning("Covariate have only 1 level. I am not estimating epsilon")
        estimate_epsilon==0
    }
    cut_points<-make_cut_points(
        onset_date=onset_date,
        household_index=household_index,
        followup_time=followup_time
    )
    ##Draw initial values for the parameters
    alpha<-runif(n=1,min=0.1,max=1)
    beta<-runif(n=1,min=0.1,max=1)
    gamma<-runif(n=1,min=0.1,max=1)
    epsilon<-rep(1,length(onset_date))
    ##Initialize the output chain
    chain<-array(data=NA,dim=c(iterations,5+2*length(levels_covariate)))
    ##Make the lookup table for hazard within the cohort.
    cohort_hazard_table<-make_cohort_hazard_table(
        onset_date=onset_date,
        followup_time=followup_time,
        constant_hazard=constant_hazard
    )
    ##Find log probability of initial parameters
    this_P<-sum(
        log(prior(
            alpha=alpha,
            beta=beta,
            gamma=gamma,
            epsilon=epsilon,
            esteps=estimate_epsilon,
            ...
            )),
        logLikelihood(
            onset=onset_date,
            cut=cut_points,
            hh=household_index,
            table=cohort_hazard_table,
            alpha=alpha,
            beta=beta,
            gamma=gamma,
            epsilon=epsilon,
            followup=followup_time,
            index=index,
            esteps=estimate_epsilon
        )
    )
    ##Initialize a vector to keep track of the acceptance rate
    accepted<-rep(0,2+length(levels_covariate))
    for(i in 1:iterations){
        ##Propose a new alpha
        proposed_alpha<-alpha*exp(rnorm(n=1,sd=delta[1]))
        proposed_P<-sum(
            log(prior(
                alpha=proposed_alpha,
                beta=beta,
                gamma=gamma,
                epsilon=epsilon,
                esteps=estimate_epsilon,
                ...
                )),
            logLikelihood(
                onset=onset_date,
                cut=cut_points,
                hh=household_index,
                table=cohort_hazard_table,
                alpha=proposed_alpha,
                beta=beta,
                gamma=gamma,
                epsilon=epsilon,
                followup=followup_time,
                index=index,
                esteps=estimate_epsilon
            )
        )
        check_probs(proposed_P,this_P)
        ##Accept or reject proposed alpha
        if(proposed_P>this_P || rbinom(n=1,size=1,prob=exp(proposed_P-this_P))==1){
            alpha<-proposed_alpha
            this_P<-proposed_P
            accepted[1]<-accepted[1]+1
        }
        ##Propose a new beta
        proposed_beta<-beta*exp(rnorm(n=1,sd=delta[2]))
        proposed_P<-sum(
            log(prior(
                alpha=alpha,
                beta=proposed_beta,
                gamma=gamma,
                epsilon=epsilon,
                esteps=estimate_epsilon,
                ...
                )),
            logLikelihood(
                onset=onset_date,
                cut=cut_points,
                hh=household_index,
                table=cohort_hazard_table,
                alpha=alpha,
                beta=proposed_beta,
                gamma=gamma,
                epsilon=epsilon,
                followup=followup_time,
                index=index,
                esteps=estimate_epsilon
            )
        )
        check_probs(proposed_P,this_P)
        ##Accept or reject proposed beta
        if(proposed_P>this_P || rbinom(n=1,size=1,prob=exp(proposed_P-this_P))==1){
            beta<-proposed_beta
            this_P<-proposed_P
            accepted[2]<-accepted[2]+1
        }
        ##Propose a new gamma
        proposed_gamma<-gamma*exp(rnorm(n=1,sd=delta[3]))
        proposed_P<-sum(
            log(prior(
                alpha=alpha,
                beta=beta,
                gamma=proposed_gamma,
                epsilon=epsilon,
                esteps=estimate_epsilon,
                ...
                )),
            logLikelihood(
                onset=onset_date,
                cut=cut_points,
                hh=household_index,
                table=cohort_hazard_table,
                alpha=alpha,
                beta=beta,
                gamma=proposed_gamma,
                epsilon=epsilon,
                followup=followup_time,
                index=index,
                esteps=estimate_epsilon
            )
        )
        check_probs(proposed_P,this_P)
        ##Accept or reject proposed gamma
        if(proposed_P>this_P || rbinom(n=1,size=1,prob=exp(proposed_P-this_P))==1){
            gamma<-proposed_gamma
            this_P<-proposed_P
            accepted[3]<-accepted[3]+1
        }
        if(estimate_epsilon==TRUE){
            k<-0
            for(X in levels_covariate[-1]){
                ##Propose new epsilon for this level of the covariate
                k<-k+1
                proposed_epsilon<-epsilon
                proposed_epsilon[which(covariate==X)]<-proposed_epsilon[which(covariate==X)]*exp(rnorm(n=1,sd=delta[4]))
                e<-which(levels_covariate==X)
                proposed_P<-sum(
                    log(prior(
                        alpha=alpha,
                        beta=beta,
                        gamma=gamma,
                        epsilon=proposed_epsilon,
                        esteps=estimate_epsilon,
                        ...
                        )),
                    logLikelihood(
                        onset=onset_date,
                        cut=cut_points,
                        hh=household_index,
                        table=cohort_hazard_table,
                        alpha=alpha,
                        beta=beta,
                        gamma=gamma,
                        epsilon=proposed_epsilon,
                        followup=followup_time,
                        index=index,
                        esteps=estimate_epsilon
                    )
                )
                check_probs(proposed_P,this_P)
                ##Accept or reject the proposed epsilon
                if(proposed_P>this_P || rbinom(n=1,size=1,prob=exp(proposed_P-this_P))==1){
                    epsilon<-proposed_epsilon
                    this_P<-proposed_P
                    accepted[3+k]<-accepted[3+k]+1
                }
            }
            chain[i,]<-c(alpha,beta,gamma,epsilon[match(x=levels_covariate[-1],table=covariate)],this_P,accepted/i)
        }else{
            chain[i,]<-c(alpha,beta,gamma,this_P,accepted/i)
        }
    }
    ##Prepare the chain for output
    if(estimate_epsilon==TRUE){
        dimnames(chain)[[2]]<-c(
            "alpha",
            "beta",
            "gamma",
            paste0("epsilon",levels_covariate[-1]),
            "logL",
            "acceptance_alpha",
            "acceptance_beta",
            "acceptance_gamma",
            paste0("acceptance_epsilon",levels_covariate[-1])
        )
    }else{
        dimnames(chain)[[2]]<-c(
            "alpha",
            "beta",
            "gamma",
            "logL",
            "acceptance_alpha",
            "acceptance_beta",
            "acceptance_gamma"
        )
    }
    class(chain)<-"chain"
    if(plot_chain){plot(chain)}
    if(estimate_epsilon==TRUE){
        cat("\nLevels of Covariate\n\n")
        print(levels_covariate)
    }
    summary(chain)
    out_list<-list(
        chain=chain,
        data=data.frame(
            onset_date=onset_date,
            household_index=household_index,
            covariate=covariate,
            epsilon=epsilon
        ),
        start_date=start_date,
        followup_time=followup_time,
        cohort_hazard_table=cohort_hazard_table,
        cut_points=cut_points,
        iterations=iterations,
        delta=delta,
        index=index,
        prior=prior
    )
    class(out_list)<-"SIRmcmc"
    return(out_list)
}

plot_SIRmcmc<-function(SIRmcmc,burnin=1,thin=1,logL=FALSE,accept=FALSE){
    plot_chain(SIRmcmc$chain,burnin=burnin,thin=thin,logL=logL,accept=accept)
}

plot_chain<-function(chain,burnin=1,thin=1,logL=FALSE,accept=FALSE){
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    plot_rows<-(ncol(chain)-1)/2+1
    plot_cols<-1
    plot_index<-seq(burnin,nrow(chain),by=thin)
    if(logL){plot_rows<-plot_rows+1}
    if(accept){plot_cols<-plot_cols+1}
    par(mfrow=c(plot_rows,plot_cols),mar=c(2,4,1,1))
    plot.default(x=plot_index,y=chain[plot_index,"alpha"],xlab="",ylab="alpha",type="l")
    plot.default(x=plot_index,y=chain[plot_index,"beta"],xlab="",ylab="beta",type="l")
    plot.default(x=plot_index,y=chain[plot_index,"gamma"],xlab="",ylab="gamma",type="l")
    plot.default(x=plot_index,y=chain[plot_index,"beta"]/chain[plot_index,"gamma"],xlab="",ylab="beta / gamma",type="l")
    if(4<=(ncol(chain)-1)/2){
        for(X in seq(4,(ncol(chain)-1)/2)){
            plot.default(x=plot_index,y=chain[plot_index,X],xlab="",ylab="epsilon",type="l")
        }
    }
    if(logL){
        plot.default(x=plot_index,y=chain[plot_index,"logL"],xlab="",ylab="log(L)",type="l")
    }
    if(accept){
        plot.default(x=plot_index,y=chain[plot_index,"acceptance_alpha"],xlab="",ylab="Accept alpha",type="l")
        plot.default(x=plot_index,y=chain[plot_index,"acceptance_beta"],xlab="",ylab="Accept beta",type="l")
        plot.default(x=plot_index,y=chain[plot_index,"acceptance_gamma"],xlab="",ylab="Accept gamma",type="l")
        if((ncol(chain)+9)/2<=ncol(chain)){
            for(X in seq((ncol(chain)+9)/2,ncol(chain))){
                plot.default(x=plot_index,y=chain[plot_index,X],xlab="",ylab="Accept epsilon",type="l")
            }
        }
    }
}

summary_SIRmcmc<-function(SIRmcmc,burnin=1,thin=1){
    summary_chain(SIRmcmc$chain,burnin=burnin,thin=thin)
}

summary_chain<-function(chain,burnin=1,thin=1){
    parameter_cols<-seq(1,(ncol(chain)-1)/2)
    accept_cols<-seq((ncol(chain)+3)/2,ncol(chain))
    print(
        t(
            apply(
                X=chain[seq(burnin,nrow(chain),thin),parameter_cols],
                MARGIN=2,
                FUN=function(x) quantile(x,probs=c(0.5,0.025,0.975)))
        )
    )
    cat("\n\n")
    print(
        chain[nrow(chain),accept_cols]
    )
}

community_attack_rate<-function(SIRmcmc,probs=c(0.5,0.025,0.975)){
    if(is(SIRmcmc, "SIRmcmc") == FALSE){stop("Object not of the SIRmcmc class.")}
    onset_date<-SIRmcmc$onset_date
    cohort_table<-SIRmcmc$cohort_hazard_table
    followup_time<-SIRmcmc$followup_time
    chain<-SIRmcmc$chain
    if(any(substr(colnames(chain),1,7)=="epsilon")){
        epsilon_cols<-which(substr(colnames(chain),1,7)=="epsilon")
        CAR<-matrix(data=NA,nrow=length(epsilon_cols)+1,ncol=length(probs))
        k<-1
        CAR[k,]<-(1-exp(-cohort_table[nrow(cohort_table)-1,"Lambda"]))*quantile(x=chain[,"alpha"],probs=probs)
        for(c in epsilon_cols){
            k<-k+1
            CAR[k,]<-(1-exp(-cohort_table[nrow(cohort_table)-1,"Lambda"]))*quantile(x=chain[,c]*chain[,"alpha"],probs=probs)
        }
        rownames(CAR)<-c("epsilon=1",colnames(chain)[epsilon_cols])
        colnames(CAR)<-probs
    }else{
        CAR<-(1-exp(-cohort_table[nrow(cohort_table)-1,"Lambda"]))*quantile(x=chain[,"alpha"],probs=probs)
        names(CAR)<-probs
    }
    return(CAR)
}

secondary_attack_rate<-function(household_size,SIRmcmc,probs=c(0.5,0.025,0.975)){
    if(is(SIRmcmc, "SIRmcmc") == FALSE){stop("Object not of the SIRmcmc class.")}
    chain<-SIRmcmc$chain
    if(any(substr(colnames(chain),1,7)=="epsilon")){
        epsilon_cols<-which(substr(colnames(chain),1,7)=="epsilon")
        SAR<-array(NA,dim=c(length(household_size),length(probs),1+length(epsilon_cols)))
        k<-1
        R_inf<-matrix(quantile(chain[,"beta"]/chain[,"gamma"],probs=probs),nrow=1,ncol=length(probs))
        SAR[,,k]<-1-exp(-matrix(1/household_size,ncol=1)%*%R_inf)
        for(c in epsilon_cols){
            k<-k+1
            R_inf<-matrix(quantile(x=chain[,c]*chain[,"beta"]/chain[,"gamma"],probs=probs),nrow=1,ncol=length(probs))
            SAR[,,k]<-1-exp(-matrix(1/household_size,ncol=1)%*%R_inf)
        }
        if(length(household_size)>1){
            dimnames(SAR)[[1]]<-household_size
        }else{
            dimnames(SAR)[[1]]<-list(household_size)
        }
        dimnames(SAR)[[2]]<-probs
        dimnames(SAR)[[3]]<-c("epsilon=1",colnames(chain)[epsilon_cols])
    }else{
        R_inf<-matrix(quantile(x=chain[,"beta"]/chain[,"gamma"],probs=probs),nrow=1,ncol=length(probs))
        SAR<-1-exp(-matrix(1/household_size,ncol=1)%*%R_inf)
        rownames(SAR)<-household_size
        colnames(SAR)<-probs
    }
    return(SAR)
}
