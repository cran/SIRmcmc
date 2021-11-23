#include <Rcpp.h> 
using namespace Rcpp;

NumericVector uniqueCpp(NumericVector x){
  NumericVector y;
  y=x[0];
  int found;
  for(NumericVector::iterator i=x.begin(); i!=x.end(); ++i){
    found=0;
    for(NumericVector::iterator j=y.begin(); j!=y.end(); ++j){
      if(*i==*j){found=1;}
    }
    if(found!=1){y.push_back(*i);}
  }
  return(y);
}

double whichHazard(double t, NumericMatrix table){
  double hazard=0;
  for(int i=0; i< table.nrow(); ++i){
    if(table(i,0)==t){
      hazard=table(i,1);
    }
  }
  return(hazard);
}


// [[Rcpp::export]]
double likelihood(
		 NumericVector onset,
		 NumericVector cut,
		 NumericMatrix table,
		 double alpha,
		 double beta,
		 double gamma,
		 NumericVector epsilon,
		 int index,
		 int esteps
		 ){
  int E;
  int n=onset.length();
  int T=cut.length();
  NumericVector uniqueEpsilon;
  if(esteps!=1){
    uniqueEpsilon=1;
    E=1;
  }
  else{
    uniqueEpsilon=uniqueCpp(epsilon);
    E=uniqueEpsilon.length();
  }
  double product=1;
  NumericVector probs;
  NumericMatrix Lambda(T,E);
  NumericMatrix LambdaPartial(T,E);
  double LambdaPartialCohort;
  double LambdaPartialHousehold;
  double a;
  double b;
  int e=0;
  double I;
  int k;
  while(e<E){
    k=1;
    while(k<T){
      a=cut(k-1);
      b=cut(k);
      I=0;
      for(NumericVector::iterator t=onset.begin(); t!=onset.end(); ++t){
	if(*t<b){
	  I +=exp(-gamma*(a-*t))-exp(-gamma*(b-*t));
	}
      }
      LambdaPartialCohort=uniqueEpsilon(e)*alpha*(whichHazard(b,table)-whichHazard(a,table));
      LambdaPartialHousehold=uniqueEpsilon(e)*beta/n/gamma*I;
      LambdaPartial(k,e)=LambdaPartialCohort+LambdaPartialHousehold;
      Lambda(k,e)=Lambda(k-1,e)+LambdaPartial(k,e);
      k += 1;
    }
    e += 1;
  }
  if(esteps<1){
    int i;
    for(NumericVector::iterator z=onset.begin(); z!=onset.end(); ++z){
      if(*z==index){
	probs.push_back(1);
      }else{
	i=1;
	while(*z!=cut(i)){
	  i += 1;
	}
	probs.push_back(exp(-Lambda(i-1,0))-exp(-Lambda(i,0)));
      }
    }
  }else{
    int i;
    int j;
    for(int z=0; z<n; ++z){
      if(onset(z)==index){
	probs.push_back(1);
      }else{
	i=1;
	while(onset(z)!=cut(i)){
	  i += 1;
	}
	j=0;
	while(epsilon(z)!=uniqueEpsilon(j)){
	  j += 1;
	}
	probs.push_back(exp(-Lambda(i-1,j))-exp(-Lambda(i,j)));
      }
    }
  }
  for(NumericVector::iterator p=probs.begin(); p!=probs.end(); ++p){
    product *= *p;
  }
  return product;
}

// [[Rcpp::export]]
double logLikelihood(
		 NumericVector onset,
		 List cut,
		 NumericVector hh,
		 NumericMatrix table,
		 double alpha,
		 double beta,
		 double gamma,
		 NumericVector epsilon,
		 int followup,
		 int index,
		 int esteps     
		     ){
  double UB=100000000000000;
  double LB=0.00000000000001;
  if((alpha<LB) || (alpha>UB) || (beta<LB) || (beta>UB) || (gamma<LB) || (gamma>UB)){
    return(-std::numeric_limits<double>::infinity());
  }
  if(esteps==1){
    for(NumericVector::iterator e=epsilon.begin(); e!=epsilon.end(); ++e){
      if((*e<LB) || (*e>UB)){
	return(-std::numeric_limits<double>::infinity());
      }
    }
  }
  long double logprob=0;
  int hhPoint=hh(0);
  NumericVector tmpOnset=NumericVector::create(onset(0));
  int tmpCut=0;
  for(int i=1; i<onset.length(); ++i){
    if(hh(i)==hhPoint){
	tmpOnset.push_back(onset(i));
    }else{
      logprob += log(likelihood(tmpOnset,cut[tmpCut],table,alpha,beta,gamma,epsilon,index,esteps));
      tmpCut += 1;
      tmpOnset=NumericVector::create();
      tmpOnset.push_back(onset(i));
      hhPoint=hh(i);
    }
  }
  logprob += log(likelihood(tmpOnset,cut[tmpCut],table,alpha,beta,gamma,epsilon,index,esteps));
  return(static_cast<double>(logprob));
}

// [[Rcpp::export]]
double unifprior(
		 double alpha,
		 double beta,
		 double gamma,
		 NumericVector epsilon,
		 int esteps,
		 double UB=1000,
		 double LB=0.01		 
		 ){
  if((alpha<LB) || (alpha>UB) || (beta<LB) || (beta>UB) || (gamma<LB) || (gamma>UB)){
    return(0);
  }
  if(esteps==1){
    for(NumericVector::iterator e=epsilon.begin(); e!=epsilon.end(); ++e){
      if((*e<LB) || (*e>UB)){
	return(0);
      }
    }
  }
  return(1);
}
