ci.blm <- function(object,C,alpha=.05,sig=4,var=NULL){

    z = qnorm(1-alpha/2)
    
    est = C%*%object@fit$par
    V = object@V

    if(!is.null(var)) V = var

    se = sqrt(t(C)%*%V%*%C)
    
    lower = est-z*se
    upper = est+z*se
		
		CI = paste(round(est,sig),", (",round(lower,sig),", ",round(upper,sig),")",sep="",coll="")

		list(est=est,se=se,lower=lower,upper=upper,CI=CI)
}

ci.lexpit <- function(object,C,alpha=.05,sig=4,baseline=TRUE,C.expit,var=NULL){

  dot.expit = function(x){exp(x)/(1+exp(x))^2}
        
        X <- model.matrix(object@formula.linear,object@data)
	int = attr(terms(object@formula.linear),"int")
	
	if(int){
		p = ncol(X)-1
		X = X[,-1]
		X = matrix(X,ncol=p)
		}
	else{
		p = ncol(X)
		}

        Z = model.matrix(object@formula.expit,object@data)
        q = ncol(Z)
	beta <- object@fit$par[1:p]
        gamma <- object@fit$par[(p+1):(p+q)]
        z = qnorm(1-alpha/2)

    if(baseline){
      if(missing(C.expit)) stop("Baseline true but linear combination C.expit not specified.")
      est = C%*%beta+expit(C.expit%*%gamma)
      C.dot = c(C,C.expit*dot.expit(C.expit%*%gamma))
    }
    else{
      est = C%*%beta
      C.dot = c(C,rep(0,q))
    }

    V = object@V
    if(!is.null(var)) V = var  

    se = sqrt(t(C.dot)%*%V%*%C.dot)
    
    lower = est-z*se
    upper = est+z*se
		
 CI = paste(round(est,sig),", (",round(lower,sig),", ",round(upper,sig),")",sep="",coll="")

 list(est=est,se=se,lower=lower,upper=upper,CI=CI)
}

setGeneric("ci",function(object,...){standardGeneric("ci")})

setMethod("ci","blm",ci.blm)
setMethod("ci","lexpit",ci.lexpit)
