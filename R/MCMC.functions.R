##Try to coerce dates to an integer format for the MCMC algorithm
MCMC.date<-function(dates,start.date=NULL){
    if(class(dates)=="character"){
        if(all(is.na(dates)==is.na(as.Date(dates)))==FALSE){
            failed.dates<-dates[is.na(dates)!=is.na(as.Date(dates))]
            failed.dates<-unique(failed.dates)
            stop(paste(
                "Coercing dates from character format to date format failed for the following:\n\n",
                failed.dates,
                "\n\nTry using as.Date() on your data."
            ))
        }else{
            dates<-as.Date(dates)
        }
    }
    if(is.null(start.date)==TRUE){
        dates<-as.numeric(dates-min(dates,na.rm=TRUE))+1
    }else{
        dates<-as.numeric(dates-start.date+1)
    }
    if(any(is.na(dates))==TRUE){
        warning("Assuming NA are non-cases")
        dates[which(is.na(dates))]<-Inf
    }
    return(dates)
}

##A lookup table for endpoints for integrating. Avoids redundant calls to unique().
make.cut.points<-function(onset.date,household.index,followup.time){
    cut.points<-list()
    length(cut.points)<-max(household.index)
    for(h in household.index){
        index<-which(household.index==h)
        temp.onset.date<-onset.date[index]
        if(all(is.finite(temp.onset.date))==FALSE){
            cut.points[[h]]<-sort(
                unique(
                    c(
                        0,
                        temp.onset.date,
                        temp.onset.date-1,
                        followup.time
                    )
                )
            )
        }else{
            cut.points[[h]]<-sort(
                unique(
                    c(
                        0,
                        temp.onset.date,
                        temp.onset.date-1
                    )
                )
            )
        }
    }
    return(cut.points)
}


##Computes cummulative hazard function within cohort. Returns a smoothed cummulative hazard and
##the product limit estimates of the cummulative hazard
make.cohort.hazard.table<-function(onset.date,followup.time,constant.hazard=FALSE){
    if(any(onset.date[which(is.finite(onset.date))]>followup.time)==TRUE){
        stop(paste("Onset date after followup ends. Rerun with followup >=",max(onset.date[which(is.finite(onset.date))])))
    }
    if(all(is.finite(onset.date))==FALSE){
        cut.points<-c(seq(0,followup.time),Inf)
    }else{
        if(all(onset.date<followup.time)==TRUE){
            stop(paste("100% attack rate before end of followup. Tried to evaluate Inf-Inf. Rerun with followup =",max(onset.date)))
        }
        cut.points<-c(seq(0,followup.time))
    }
    K<-length(cut.points)
    cases<-table(
        factor(
            x=onset.date,
            levels=cut.points
        )
    )
    survival<-rep(1,K)
    for(k in 2:K){
        survival[k]<-survival[k-1]*(1-cases[k]/sum(cases[k:K]))
    }
    Lambda<--log(survival)
    lambda<-diff(Lambda)
    if(length(K)<=7){
        if(class(try(loess(lambda[1:(K-2)]~cut.points[1:(K-2)])))!="try-error"){
            smooth.lambda<-loess(lambda[1:(K-2)]~cut.points[1:(K-2)])$fitted
            if(any(smooth.lambda<0)){
                smooth.lambda<-loess(lambda[1:(K-2)]~cut.points[1:(K-2)],degree=0)$fitted
                warning("Using degree 0 loess")
            }
            if(any(smooth.lambda<0)){
                smooth.lambda<-lambda
                warning("Using empiric survival")
            }
        }else{
            smooth.lambda<-lambda
            warning("LOESS failed. Using empiric suvival")
        }
    }else{
        smooth.lambda<-lambda
        warning("Too few time points. Using empiric survival.")
    }
    smooth.Lambda<-cumsum(smooth.lambda)
    if(constant.hazard==FALSE){
        return(
            cbind(
                t=cut.points,
                Lambda=c(0,smooth.Lambda,Inf),
                empiric=Lambda
            )
        )
    }else{
        return(
            cbind(
                t=cut.points,
                Lambda=cut.points,
                empiric=Lambda
            )
        )
    }
}

##Sanity checks on log probabilities
check.probs<-function(p,q){
    if(is.na(p)==TRUE | is.na(q)){stop("Probability NA")}
    if(is.finite(p)==FALSE & is.finite(q)==FALSE){stop("Tried to Evaluate -Inf+Inf")}
    if(is.numeric(p)==FALSE | is.numeric(q)==FALSE){stop("Probability not a numeric variable")}
    if(is.nan(p)==TRUE | is.nan(q)==TRUE){stop(paste("Probability not a number"))}
}

##Find the index that sorts households, smallest households first.
sort.order.households<-function(household.index){
    n<-length(household.index)
    out.index<-rep(NA,length(n))
    hh.table<-sort(table(household.index))
    i<-1
    for(h in as.numeric(names(hh.table))){
        temp.index<-which(household.index==h)
        m<-length(temp.index)
        out.index[i:(i+m-1)]<-temp.index
        i<-i+m
    }
    return(out.index)
}

##Renumber the households in increasing order
renumber.households<-function(household.index){
    n<-length(household.index)
    out.index<-rep(NA,length(n))
    new.index<-1
    out.index[1]<-new.index
    for(i in 2:n){
        if(household.index[i]!=household.index[i-1]){
            new.index<-new.index+1
        }
        out.index[i]<-new.index
    }
    return(out.index)
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
household.transmission<-function(
                                 onset.date,
                                 household.index,
                                 covariate=NULL,
                                 followup.time,
                                 iterations,
                                 delta=c(0.1,0.3,0.4,0.1),
                                 plot.chain=TRUE,
                                 index=1,
                                 start.date=NULL,
                                 prior=unifprior,
                                 constant.hazard=FALSE,
                                 ...
                                 ){
    ##Preprares the data for the algorithm: households together, small households first, and an increasing index.
    household.index<-as.integer(as.factor(household.index))
    onset.date<-MCMC.date(onset.date,start.date)
    sort.order<-sort.order.households(household.index)
    onset.date<-onset.date[sort.order]
    household.index<-household.index[sort.order]
    household.index<-renumber.households(household.index)
    ##Decide whether to estimate epsilon, and structure the covariate data for the algorithm.
    if(is.null(covariate)){
        covariate<-rep(1,length(onset.date))
        estimate.epsilon<-0
    }else{
        estimate.epsilon<-1
        covariate<-covariate[sort.order]
    }
    levels.covariate<-sort(unique(covariate))
    if(length(levels.covariate)==1 & estimate.epsilon==1){
        warning("Covariate have only 1 level. I am not estimating epsilon")
        estimate.epsilon==0
    }
    cut.points<-make.cut.points(
        onset.date=onset.date,
        household.index=household.index,
        followup.time=followup.time
    )
    ##Draw initial values for the parameters
    alpha<-runif(n=1,min=0.1,max=1)
    beta<-runif(n=1,min=0.1,max=1)
    gamma<-runif(n=1,min=0.1,max=1)
    epsilon<-rep(1,length(onset.date))
    ##Initialize the output chain
    chain<-array(data=NA,dim=c(iterations,5+2*length(levels.covariate)))
    ##Make the lookup table for hazard within the cohort.
    cohort.hazard.table<-make.cohort.hazard.table(
        onset.date=onset.date,
        followup.time=followup.time,
        constant.hazard=constant.hazard
    )
    ##Find log probability of initial parameters
    this.P<-sum(
        log(prior(
            alpha=alpha,
            beta=beta,
            gamma=gamma,
            epsilon=epsilon,
            esteps=estimate.epsilon,
            ...
            )),
        logLikelihood(
            onset=onset.date,
            cut=cut.points,
            hh=household.index,
            table=cohort.hazard.table,
            alpha=alpha,
            beta=beta,
            gamma=gamma,
            epsilon=epsilon,
            followup=followup.time,
            index=index,
            esteps=estimate.epsilon
        )
    )
    ##Initialize a vector to keep track of the acceptance rate
    accepted<-rep(0,2+length(levels.covariate))
    for(i in 1:iterations){
        ##Propose a new alpha
        proposed.alpha<-alpha*exp(rnorm(n=1,sd=delta[1]))
        proposed.P<-sum(
            log(prior(
                alpha=proposed.alpha,
                beta=beta,
                gamma=gamma,
                epsilon=epsilon,
                esteps=estimate.epsilon,
                ...
                )),
            logLikelihood(
                onset=onset.date,
                cut=cut.points,
                hh=household.index,
                table=cohort.hazard.table,
                alpha=proposed.alpha,
                beta=beta,
                gamma=gamma,
                epsilon=epsilon,
                followup=followup.time,
                index=index,
                esteps=estimate.epsilon
            )
        )
        check.probs(proposed.P,this.P)
        ##Accept or reject proposed alpha
        if(proposed.P>this.P || rbinom(n=1,size=1,prob=exp(proposed.P-this.P))==1){
            alpha<-proposed.alpha
            this.P<-proposed.P
            accepted[1]<-accepted[1]+1
        }
        ##Propose a new beta
        proposed.beta<-beta*exp(rnorm(n=1,sd=delta[2]))
        proposed.P<-sum(
            log(prior(
                alpha=alpha,
                beta=proposed.beta,
                gamma=gamma,
                epsilon=epsilon,
                esteps=estimate.epsilon,
                ...
                )),
            logLikelihood(
                onset=onset.date,
                cut=cut.points,
                hh=household.index,
                table=cohort.hazard.table,
                alpha=alpha,
                beta=proposed.beta,
                gamma=gamma,
                epsilon=epsilon,
                followup=followup.time,
                index=index,
                esteps=estimate.epsilon
            )
        )
        check.probs(proposed.P,this.P)
        ##Accept or reject proposed beta
        if(proposed.P>this.P || rbinom(n=1,size=1,prob=exp(proposed.P-this.P))==1){
            beta<-proposed.beta
            this.P<-proposed.P
            accepted[2]<-accepted[2]+1
        }
        ##Propose a new gamma
        proposed.gamma<-gamma*exp(rnorm(n=1,sd=delta[3]))
        proposed.P<-sum(
            log(prior(
                alpha=alpha,
                beta=beta,
                gamma=proposed.gamma,
                epsilon=epsilon,
                esteps=estimate.epsilon,
                ...
                )),
            logLikelihood(
                onset=onset.date,
                cut=cut.points,
                hh=household.index,
                table=cohort.hazard.table,
                alpha=alpha,
                beta=beta,
                gamma=proposed.gamma,
                epsilon=epsilon,
                followup=followup.time,
                index=index,
                esteps=estimate.epsilon
            )
        )
        check.probs(proposed.P,this.P)
        ##Accept or reject proposed gamma
        if(proposed.P>this.P || rbinom(n=1,size=1,prob=exp(proposed.P-this.P))==1){
            gamma<-proposed.gamma
            this.P<-proposed.P
            accepted[3]<-accepted[3]+1
        }
        if(estimate.epsilon==TRUE){
            k<-0
            for(X in levels.covariate[-1]){
                ##Propose new epsilon for this level of the covariate
                k<-k+1
                proposed.epsilon<-epsilon
                proposed.epsilon[which(covariate==X)]<-proposed.epsilon[which(covariate==X)]*exp(rnorm(n=1,sd=delta[4]))
                e<-which(levels.covariate==X)
                proposed.P<-sum(
                    log(prior(
                        alpha=alpha,
                        beta=beta,
                        gamma=gamma,
                        epsilon=proposed.epsilon,
                        esteps=estimate.epsilon,
                        ...
                        )),
                    logLikelihood(
                        onset=onset.date,
                        cut=cut.points,
                        hh=household.index,
                        table=cohort.hazard.table,
                        alpha=alpha,
                        beta=beta,
                        gamma=gamma,
                        epsilon=proposed.epsilon,
                        followup=followup.time,
                        index=index,
                        esteps=estimate.epsilon
                    )
                )
                check.probs(proposed.P,this.P)
                ##Accept or reject the proposed epsilon
                if(proposed.P>this.P || rbinom(n=1,size=1,prob=exp(proposed.P-this.P))==1){
                    epsilon<-proposed.epsilon
                    this.P<-proposed.P
                    accepted[3+k]<-accepted[3+k]+1
                }
            }
            chain[i,]<-c(alpha,beta,gamma,epsilon[match(x=levels.covariate[-1],table=covariate)],this.P,accepted/i)
        }else{
            chain[i,]<-c(alpha,beta,gamma,this.P,accepted/i)
        }
    }
    ##Prepare the chain for output
    if(estimate.epsilon==TRUE){
        dimnames(chain)[[2]]<-c(
            "alpha",
            "beta",
            "gamma",
            paste0("epsilon",levels.covariate[-1]),
            "logL",
            "acceptance.alpha",
            "acceptance.beta",
            "acceptance.gamma",
            paste0("acceptance.epsilon",levels.covariate[-1])
        )
    }else{
        dimnames(chain)[[2]]<-c(
            "alpha",
            "beta",
            "gamma",
            "logL",
            "acceptance.alpha",
            "acceptance.beta",
            "acceptance.gamma"
        )
    }
    class(chain)<-"chain"
    if(plot.chain){plot(chain)}
    if(estimate.epsilon==TRUE){
        cat("\nLevels of Covariate\n\n")
        print(levels.covariate)
    }
    summary(chain)
    out.list<-list(
        chain=chain,
        data=data.frame(
            onset.date=onset.date,
            household.index=household.index,
            covariate=covariate,
            epsilon=epsilon
        ),
        start.date=start.date,
        followup.time=followup.time,
        cohort.hazard.table=cohort.hazard.table,
        cut.points=cut.points,
        iterations=iterations,
        delta=delta,
        index=index,
        prior=prior
    )
    class(out.list)<-"SIRmcmc"
    return(out.list)
}

plot.SIRmcmc<-function(SIRmcmc,burnin=1,thin=1,logL=FALSE,accept=FALSE){
    plot.chain(SIRmcmc$chain,burnin=burnin,thin=thin,logL=logL,accept=accept)
}

plot.chain<-function(chain,burnin=1,thin=1,logL=FALSE,accept=FALSE){
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    plot.rows<-(ncol(chain)-1)/2+1
    plot.cols<-1
    plot.index<-seq(burnin,nrow(chain),by=thin)
    if(logL){plot.rows<-plot.rows+1}
    if(accept){plot.cols<-plot.cols+1}
    par(mfrow=c(plot.rows,plot.cols),mar=c(2,4,1,1))
    plot.default(x=plot.index,y=chain[plot.index,"alpha"],xlab="",ylab="\U03B1",type="l")
    plot.default(x=plot.index,y=chain[plot.index,"beta"],xlab="",ylab="\U03B2",type="l")
    plot.default(x=plot.index,y=chain[plot.index,"gamma"],xlab="",ylab="\U03B3",type="l")
    plot.default(x=plot.index,y=chain[plot.index,"beta"]/chain[plot.index,"gamma"],xlab="",ylab="\U03B2 / \U03B3",type="l")
    if(4<=(ncol(chain)-1)/2){
        for(X in seq(4,(ncol(chain)-1)/2)){
            plot.default(x=plot.index,y=chain[plot.index,X],xlab="",ylab="\U03B5",type="l")
        }
    }
    if(logL){
        plot.default(x=plot.index,y=chain[plot.index,"logL"],xlab="",ylab="log(L)",type="l")
    }
    if(accept){
        plot.default(x=plot.index,y=chain[plot.index,"acceptance.alpha"],xlab="",ylab="Accept \U03B1",type="l")
        plot.default(x=plot.index,y=chain[plot.index,"acceptance.beta"],xlab="",ylab="Accept \U03B2",type="l")
        plot.default(x=plot.index,y=chain[plot.index,"acceptance.gamma"],xlab="",ylab="Accept\U03B3",type="l")
        if((ncol(chain)+9)/2<=ncol(chain)){
            for(X in seq((ncol(chain)+9)/2,ncol(chain))){
                plot.default(x=plot.index,y=chain[plot.index,X],xlab="",ylab="Accept \U03B5",type="l")
            }
        }
    }
}

summary.SIRmcmc<-function(SIRmcmc,burnin=1,thin=1){
    summary.chain(SIRmcmc$chain,burnin=burnin,thin=thin)
}

summary.chain<-function(chain,burnin=1,thin=1){
    parameter.cols<-seq(1,(ncol(chain)-1)/2)
    accept.cols<-seq((ncol(chain)+3)/2,ncol(chain))
    print(
        t(
            apply(
                X=chain[seq(burnin,nrow(chain),thin),parameter.cols],
                MARGIN=2,
                FUN=function(x) quantile(x,probs=c(0.5,0.025,0.975)))
        )
    )
    cat("\n\n")
    print(
        chain[nrow(chain),accept.cols]
    )
}

community.attack.rate<-function(SIRmcmc,probs=c(0.5,0.025,0.975)){
    if(class(SIRmcmc) != "SIRmcmc"){stop("Object not of the SIRmcmc class.")}
    onset.date<-SIRmcmc$onset.date
    cohort.table<-SIRmcmc$cohort.hazard.table
    followup.time<-SIRmcmc$followup.time
    chain<-SIRmcmc$chain
    if(any(substr(colnames(chain),1,7)=="epsilon")){
        epsilon.cols<-which(substr(colnames(chain),1,7)=="epsilon")
        CAR<-matrix(data=NA,nrow=length(epsilon.cols)+1,ncol=length(probs))
        k<-1
        CAR[k,]<-(1-exp(-cohort.table[nrow(cohort.table)-1,"Lambda"]))*quantile(x=chain[,"alpha"],probs=probs)
        for(c in epsilon.cols){
            k<-k+1
            CAR[k,]<-(1-exp(-cohort.table[nrow(cohort.table)-1,"Lambda"]))*quantile(x=chain[,c]*chain[,"alpha"],probs=probs)
        }
        rownames(CAR)<-c("epsilon=1",colnames(chain)[epsilon.cols])
        colnames(CAR)<-probs
    }else{
        CAR<-(1-exp(-cohort.table[nrow(cohort.table)-1,"Lambda"]))*quantile(x=chain[,"alpha"],probs=probs)
        names(CAR)<-probs
    }
    return(CAR)
}

secondary.attack.rate<-function(household.size,SIRmcmc,probs=c(0.5,0.025,0.975)){
    if(class(SIRmcmc) != "SIRmcmc"){stop("Object not of the SIRmcmc class.")}
    chain<-SIRmcmc$chain
    if(any(substr(colnames(chain),1,7)=="epsilon")){
        epsilon.cols<-which(substr(colnames(chain),1,7)=="epsilon")
        SAR<-array(NA,dim=c(length(household.size),length(probs),1+length(epsilon.cols)))
        k<-1
        R.inf<-matrix(quantile(chain[,"beta"]/chain[,"gamma"],probs=probs),nrow=1,ncol=length(probs))
        SAR[,,k]<-1-exp(-matrix(1/household.size,ncol=1)%*%R.inf)
        for(c in epsilon.cols){
            k<-k+1
            R.inf<-matrix(quantile(x=chain[,c]*chain[,"beta"]/chain[,"gamma"],probs=probs),nrow=1,ncol=length(probs))
            SAR[,,k]<-1-exp(-matrix(1/household.size,ncol=1)%*%R.inf)
        }
        if(length(household.size)>1){
            dimnames(SAR)[[1]]<-household.size
        }else{
            dimnames(SAR)[[1]]<-list(household.size)
        }
        dimnames(SAR)[[2]]<-probs
        dimnames(SAR)[[3]]<-c("epsilon=1",colnames(chain)[epsilon.cols])
    }else{
        R.inf<-matrix(quantile(x=chain[,"beta"]/chain[,"gamma"],probs=probs),nrow=1,ncol=length(probs))
        SAR<-1-exp(-matrix(1/household.size,ncol=1)%*%R.inf)
        rownames(SAR)<-household.size
        colnames(SAR)<-probs
    }
    return(SAR)
}
