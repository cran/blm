predict.blm <- function(object,...){

  X <- model.matrix(object@formula,object@data)
  Y <- model.frame(object@formula,object@data)[,1]
  X.scaled <- apply(X,2,scale)

  #PREDICTED VALUES
  X%*%object@fit$par
}

predict.lexpit <- function(object,...){

  X <- model.matrix(object@formula.linear,object@data)
  Y <- model.frame(object@formula.linear,object@data)[,1]
  Z <- model.matrix(object@formula.expit,object@data)

  if(attr(terms(object@formula.linear),"intercept")==1) X = X[,-1]
  if(!is.matrix(X)) X = matrix(X,ncol=1)
 
  #PREDICTED VALUES
  beta = summary(object)$est.linear
  gamma = summary(object)$est.expit
  
  X%*%beta+expit(Z%*%gamma)
}

setMethod("predict",signature(object="blm"),predict.blm)
setMethod("predict",signature(object="lexpit"),predict.lexpit)


resid.blm <- function(object,...){

  Y <- model.frame(object@formula,object@data)[,1]
  Y-predict(object)
}

resid.lexpit <- function(object,...){

  Y <- model.frame(object@formula.linear,object@data)[,1]
  Y-predict(object)
}

setMethod("resid",signature(object="blm"),resid.blm)
setMethod("resid",signature(object="lexpit"),resid.lexpit)


ci.for.blm <- function(object,C,alpha=.05,subsample=FALSE,sig=4,...){

    z = qnorm(1-alpha/2)
    est = C%*%object@fit$par

    if(!subsample){
      
    V = object@V
    se = sqrt(t(C)%*%V%*%C)
    
    lower = est-z*se
    upper = est+z*se
		
  }
  else{
    B = subresample(object,...)
    E = apply(B,1,function(x){sum(C*x)})
    lower = quantile(E,alpha/2)
    upper = quantile(E,1-alpha/2)
    se = sd(E)
   }

   CI = paste(round(est,sig),", (",round(lower,sig),", ",round(upper,sig),")",sep="",coll="")   
   list(est=est,se=se,lower=lower,upper=upper,CI=CI)
}

blm.ci <- function(object,C,alpha=.05,
                   sig=4,
                   method=c("hessian","subsample","logit"),
                   ...){

   method <- match.arg(method)

   if(is.constrained(object)&method=="hessian"){
     warning("Large-sample normal CI might be inaccurate when active constraints are present.")
   }
   
   if(method!="logit"){
      ci.blm(object,C,alpha=alpha,sig=sig,...)
    }
    else{
      ci.logit.transform(object,C,alpha=alpha,sig=sig,...)
    }
}

ci.for.lexpit <- function(object,C,alpha=.05,subsample=FALSE,baseline=TRUE,sig=4,...){

if(!is.list(C)){
  stop("C must be a list with linear and expit arguments for lexpit model.")
}

C.linear = C$linear
C.expit = C$expit

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
        z = qnorm(1-alpha/2)
  
 get.est <- function(par,C,C.expit){

        beta <- par[1:p]
        gamma <- par[(p+1):(p+q)]
   
    if(baseline){
      if(missing(C.expit)) stop("Baseline true but linear combination C.expit not specified.")
      est = C%*%beta+expit(C.expit%*%gamma)
      C.dot = c(C,C.expit*dot.expit(C.expit%*%gamma))
    }
    else{
      est = C%*%beta
      C.dot = c(C,rep(0,q))
    }

  return(list(est=est,C.dot=C.dot))
}
    est.obj <- get.est(object@fit$par,C.linear,C.expit)
    est = est.obj$est

 if(!subsample){

    V = object@V
    
    se = sqrt(t(est.obj$C.dot)%*%V%*%est.obj$C.dot)
    
    lower = est-z*se
    upper = est+z*se
}
 else{
    B = subresample(object,...)
    E = apply(B,1,function(x){get.est(x,C.linear,C.expit)$est})
    lower = quantile(E,alpha/2)
    upper = quantile(E,1-alpha/2)
    se = sd(E)
   }
		
 CI = paste(round(est,sig),", (",round(lower,sig),", ",round(upper,sig),")",sep="",coll="")
 list(est=est,se=se,lower=lower,upper=upper,CI=CI)
}

setGeneric("ci.blm",function(object,...){standardGeneric("ci.blm")})
setMethod("ci.blm","blm",ci.for.blm)
setMethod("ci.blm","lexpit",ci.for.lexpit)


ci.logit.blm <- function(object,C,alpha=.05,coef=TRUE,average,sig=4){

g1 = g0 = expit
g1.dot = g0.dot = function(x){exp(x)/(1+exp(x))^2}

C1.index = which(!average)

X <- model.matrix(object@formula,object@data)
fit.glm <- glm(object@formula,object@data,family="binomial")

beta <- coef(fit.glm)
V <- vcov(fit.glm)
z = qnorm(1 - alpha/2)

X1 <- X

for(i in C1.index){
 X1[,i] <- C[i]
}

covariate.class <- apply(X1,1,function(x){paste(x,sep="",collapse="")})
covariate.patterns <- unique(covariate.class) #IN SAME ORDER AS unique(X)
p = sapply(covariate.patterns,function(x){mean(covariate.class==x)})

X1.bar <- unique(X1)
X0.bar <- X1.bar

if(!coef){

pi1 = g1(X1.bar%*%beta)
logit.beta = sum(p*pi1)
pi1.dot <- g1.dot(pi1)

X.dot <- X1.bar*as.numeric(pi1.dot)

v <- apply(X.dot,1,function(x){t(x)%*%V%*%x})
se <- sqrt(sum(p^2*v))

}
else{

for(i in C1.index){
 X0.bar[,i] <- 0
}

pi1 = g1(X1.bar%*%beta)
pi0 = g1(X0.bar%*%beta)

logit.beta = sum(p*(pi1-pi0))

pi1.dot <- g1.dot(pi1)
pi0.dot <- g0.dot(pi0)

X.dot <- X1.bar*as.numeric(pi1.dot)-X0.bar*as.numeric(pi0.dot)

v <- apply(X.dot,1,function(x){t(x)%*%V%*%x})
se <- sqrt(sum(p^2*v))

}

logit.beta = C%*%coef(object)
lower = logit.beta-z*se
upper = logit.beta+z*se
CI = paste(round(logit.beta,sig),", (",round(lower,sig),", ",round(upper,sig),")",sep="",coll="")

list(
     est = logit.beta,
     lower = lower,
     upper = upper,
     CI = CI
     )
}


ci.logit.lexpit <- function(object,C,alpha=.05,sig=4,coef=TRUE,average){

if(!is.list(C)){
  stop("C and average must be a list with linear and expit arguments for lexpit model.")
}

C1 = c(C$expit[1],C$linear,C$expit[-1])
average = c(average$expit[1],average$linear,average$expit[-1])
C1.index = which(!average)
  
g1 = g0 = expit
g1.dot = g0.dot = function(x){exp(x)/(1+exp(x))^2}

f.expit <- paste("~.+",paste(object@formula.expit)[3],collapse="")
f <- update(object@formula.linear,f.expit)
fit.glm <- glm(f,object@data,family="binomial")
X <- model.matrix(f,object@data)
  
beta <- coef(fit.glm)
V <- vcov(fit.glm)
z = qnorm(1 - alpha/2)

X1 <- X

for(i in C1.index){
 X1[,i] <- C1[i]
}

covariate.class <- apply(X1,1,function(x){paste(x,sep="",collapse="")})
covariate.patterns <- unique(covariate.class) #IN SAME ORDER AS unique(X)
p = sapply(covariate.patterns,function(x){mean(covariate.class==x)})

X1.bar <- unique(X1)
X0.bar <- X1.bar

p.linear = length(C$linear)
p.expit = length(C$expit)

par.linear = coef(object)[1:p.linear]
par.expit = coef(object)[(p.linear+1):(p.linear+p.expit)]

if(!coef){

pi1 = g1(X1.bar%*%beta)

logit.beta = sum(p*pi1)
pi1.dot <- g1.dot(pi1)

X.dot <- X1.bar*as.numeric(pi1.dot)

v <- apply(X.dot,1,function(x){t(x)%*%V%*%x})
se <- sqrt(sum(p^2*v))

logit.beta = sum(C$linear*par.linear)+expit(sum(C$expit*par.expit))

}
else{

  if(any(C$expit==1)){ #NO AVERAGING FOR EXPIT TERMS
    i <- which(C1!=0)
    logit.beta = sum(C$expit*par.expit)
    se <- sqrt(V[i,i])    
  }
  else{
    
    for(i in which(C1==0&average)){
      X0.bar[,i] <- 0
     }

 pi1 = g1(X1.bar%*%beta)
 pi0 = g1(X0.bar%*%beta)

 logit.beta = sum(p*(pi1-pi0))

 pi1.dot <- g1.dot(pi1)
 pi0.dot <- g0.dot(pi0)

 X.dot <- X1.bar*as.numeric(pi1.dot)-X0.bar*as.numeric(pi0.dot)

 v <- apply(X.dot,1,function(x){t(x)%*%V%*%x})
 se <- sqrt(sum(p^2*v))

 logit.beta = sum(C$linear*par.linear)
  }
}

lower = logit.beta-z*se
upper = logit.beta+z*se
CI = paste(round(logit.beta,sig),", (",round(lower,sig),", ",round(upper,sig),")",sep="",coll="")

list(
     est = logit.beta,
     lower = lower,
     upper = upper,
     CI = CI
     )
}

setGeneric("ci.logit.transform",function(object,C,alpha=0.05,...){standardGeneric("ci.logit.transform")})
setMethod("ci.logit.transform","blm",ci.logit.blm)
setMethod("ci.logit.transform","lexpit",ci.logit.lexpit)


setGeneric("ci",function(object,C,alpha=0.05,sig=4,method=c("hessian","subsample","logit"),...){standardGeneric("ci")})
setMethod("ci","blm",blm.ci)
setMethod("ci","lexpit",blm.ci)

setMethod("vcov","blm",function(object,...){object@V})
setMethod("vcov","lexpit",function(object,...){object@V})

setMethod("summary",signature(object="blm"),function(object){

	X <- model.matrix(object@formula,object@data)
        Y <- model.frame(object@formula,object@data)[,1]
        N <- nrow(object@data)
        loglik.null <- glm(Y~1,family=binomial)
        loglik.null <- -(summary(loglik.null)$aic-2)/2
        loglik <- -object@fit$value
	pre <- X%*%object@fit$par        
	beta = object@fit$par
	p = nrow(beta)
	rownames(beta) = colnames(X)
	
	list(
				est=beta,
				gradient=object@f.score(beta),
				feasible=all(pre>0&pre<1),
                                active = object@active.constraints$active.constraints,
				convergence=object@fit$con,
				message=object@fit$mess,
				loglik= loglik,
                                loglik.null=loglik.null,
				df = nrow(X)-p,
				AIC = 2*p +2*object@fit$value,
				null.deviance = 2*object@f.loglik(c(beta[1],rep(0,p-1))),
                                R2.efron = efron(Y,pre),
                                R2.mcfadden = mcfadden(loglik,loglik.null),
                                R2.mcfadden.adj = mcfadden.adj(loglik,loglik.null,p),
                                R2.coxsnell = cox.snell(loglik,loglik.null,N),
                                R2.coxsnell.adj = cox.snell.adj(loglik,loglik.null,N), 
                                seconds.to.run=object@run.time
	)

})

setMethod("summary",signature(object="lexpit"),function(object){

	X <- model.matrix(object@formula.linear,object@data)
	int = attr(terms(object@formula.linear),"int")
        beta.names <- colnames(object@par.start$names)

	if(int){
		p = ncol(X)-1
		X = X[,-1]
		X = matrix(X,ncol=p)
		}
	else{
		p = ncol(X)
		}

	
	Z <- model.matrix(object@formula.expit,object@data)
	q <- ncol(Z)
	
	beta <- object@fit$par[1:p]
	gamma <- object@fit$par[(p+1):(p+q)]

        Y <- model.frame(object@formula.linear,object@data)[,1]
        N <- nrow(object@data)
        loglik.null <- glm(Y~1,family=binomial)
        loglik.null <- -(summary(loglik.null)$aic-2)/2
        loglik <- -object@fit$value
	pre <- X%*%beta+expit(Z%*%gamma)
	
	names(beta) <- beta.names
	names(gamma) <- colnames(Z)
	
	list(           	est.linear=beta,
				est.expit=gamma,				
				baseline.risk=expit(gamma[1]),
				OR=exp(gamma[-1]),
				gradient=object@f.score(object@fit$par),
				feasible=all(pre>0&pre<1),
                                active = object@active.constraints$active.constraints,
				convergence=object@fit$con,
				message=object@fit$mess,
				loglik=loglik,
                                loglik.null=loglik.null,
				df = nrow(X)-(p+q),
				AIC = 2*(p+q) +2*object@fit$value,
				null.deviance = 2*object@f.loglik(c(rep(0,p),gamma[1],rep(0,q-1))),
                                R2.efron = efron(Y,pre),
                                R2.mcfadden = mcfadden(loglik,loglik.null),
                                R2.mcfadden.adj = mcfadden.adj(loglik,loglik.null,p+q),
                                R2.coxsnell = cox.snell(loglik,loglik.null,N),
                                R2.coxsnell.adj = cox.snell.adj(loglik,loglik.null,N), 
                                seconds.to.run = object@run.time
	)

})

boot.blm <- function(object,n.boot=25){

		mle = object@fit$par
		p = length(mle)
		boot.estimates <- matrix(0,n.boot,p)
		N <- nrow(object@data)

                i = 0      
		while(i < n.boot){

                  sample <- sample(1:N,replace=TRUE)
		  data <- object@data[sample,]
                                
		  one.step.fit <-  tryCatch(blm(object@formula,data,	   	                                                   ineq = object@ineq,
                                par.init=mle,
                                control.outer=list(trace=FALSE)),
                                error=function(e){NA})
                  
                  if(class(one.step.fit)!="logical"){
                    i = i + 1
                    boot.estimates[i,] <- one.step.fit@fit$par
                   }
		}
                                
   boot.estimates
}

boot.lexpit <- function(object,n.boot=25){

				mle = object@fit$par
                                fit = summary(object)
                                gamma = fit$est.expit
                                beta = fit$est.linear
				p = length(mle)
				boot.estimates <- matrix(0,n.boot,p)
				N <- nrow(object@data)

                i = 0      
		while(i < n.boot){
		  sample <- sample(1:N,replace=TRUE)
		  data <- object@data[sample,]
                                
		  one.step.fit <-  tryCatch(lexpit(object@formula.linear,object@formula.expit,data,	   	                          ineq = object@ineq,
                                       par.init=list(linear=beta,expit=gamma),
                                       control.outer=list(trace=FALSE)),
                                            error=function(e){NA})
                  
                  if(class(one.step.fit)!="logical"){
                    i = i + 1
                    boot.estimates[i,] <- one.step.fit@fit$par
                   }
		}
				
  boot.estimates
}

subboot.blm <- function(object,m=NULL,n.boot=25){

		mle = object@fit$par
		p = length(mle)
		boot.estimates <- matrix(0,n.boot,p)
		N <- nrow(object@data)
                w <- object@weights
                weighted <- ifelse(all(w==1),FALSE,TRUE)
                
                i = 0      
		while(i < n.boot){
                  
                  if(is.null(m)){
                    min.per <- ifelse(10*p/N>1,1,10*p/N)                    
                    sample.size <- ceiling(N^(runif(1,log(N*min.per)/log(N),1)))
                  }
                  else{
                    sample.size = m
                  }

                 
                  replace = ifelse(sample.size==N,TRUE,FALSE)
                  sample <- sample(1:N,size=sample.size,replace=replace)
		  data <- object@data[sample,]

                  ws <- object@weights[sample]
                  if(!weighted) ws <- NULL
                  
		  one.step.fit <-  tryCatch(blm(object@formula,data,	   	                                                    ineq = object@ineq,
                                par.init=mle,weights=ws,
                                control.outer=list(trace=FALSE)),
                                error=function(e){NA})
                  
                  if(class(one.step.fit)!="logical"){
                    i = i + 1
                    boot.estimates[i,] <- one.step.fit@fit$par
                   }
		}
                                
   boot.estimates
}

subboot.lexpit <- function(object,m=NULL,n.boot=25){

				mle = object@fit$par
                                fit = summary(object)
                                gamma = fit$est.expit
                                beta = fit$est.linear
				p = length(mle)
				boot.estimates <- matrix(0,n.boot,p)
				N <- nrow(object@data)

                w <- object@weights
                weighted <- ifelse(all(w==1),FALSE,TRUE)

                i = 0      
		while(i < n.boot){
                  
                  
                  if(is.null(m)){
                    min.per <- ifelse(10*p/N>1,1,10*p/N)                    
                    sample.size <- ceiling(N^(runif(1,log(N*min.per)/log(N),1)))
                  }
                  else{
                    sample.size = m
                  }

                  replace = ifelse(sample.size==N,TRUE,FALSE)
                  sample <- sample(1:N,size=sample.size,replace=replace)
                  ws <- object@weights[sample]

                  if(!weighted) ws <- NULL
                  
		  data <- object@data[sample,]
                                
		  one.step.fit <-  tryCatch(lexpit(object@formula.linear,object@formula.expit,
                                       data,ineq = object@ineq,weights=ws,
                                       par.init=list(linear=beta,expit=gamma),
                                       control.outer=list(trace=FALSE)),
                                            error=function(e){NA})
                  
                  if(class(one.step.fit)!="logical"){
                    i = i + 1
                    boot.estimates[i,] <- one.step.fit@fit$par
                   }
		}
				
  boot.estimates
}



setGeneric("subresample",function(object,m=NULL,n.boot=25){standardGeneric("subresample")})
setMethod("subresample","blm",subboot.blm)
setMethod("subresample","lexpit",subboot.lexpit)

print.like.glm <- function(x,digits=4){

  ## MODIFIED FROM print.glm
     beta = coef(x)
     p = nrow(beta)
     df = nrow(x@data)-p
     
     if(class(x)[1]=="blm"){
     cat("\nCall:  ", paste(deparse(x@formula), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
   }
     else{
     cat("\nLinear Call:  ", paste(deparse(x@formula.linear), sep = "\n", collapse = "\n"),
         sep = "")
          cat("\nExpit Call:  ", paste(deparse(x@formula.expit), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
     }
        cat("Coefficients")
        cat(":\n")

        T <- beta/sqrt(diag(vcov(x)))
        V <- sqrt(diag(vcov(x)))
        P <- mapply(pt,q=abs(T),MoreArgs=list(df=df))
        P <- 2*(1-P)
        coeftable <- cbind(beta,T,V,P)
        coeftable <- format(coeftable,digits=digits,sci=TRUE)
        colnames(coeftable) <- c("estimates","t-value","std. err","p-value")
                       
        print.default(coeftable, 
        print.gap = 2, quote = FALSE)

    cat("\nDegrees of Freedom:",df,"\n")
    cat("Run time (sec):",x@run.time,"\n")
    LL <- -x@f.loglik(beta)
    AIC <- 2*length(beta)-2*LL
    cat("LogLik:\t   ", format(signif(LL, 
        digits)), "\tAIC:", format(signif(AIC, digits)))
                 
    invisible(x)
}


setMethod("print","blm",print.like.glm)
setMethod("print","lexpit",print.like.glm)

setMethod("show",signature(object="blm"),function(object){

               print(object)
  
})

setMethod("show",signature(object="lexpit"),function(object){

               print(object)
  
})

gof.blm <- function(object,groups=10){

  Y <- model.frame(object@formula,object@data)[,1]

  #PREDICTED VALUES
  prediction <- predict(object)
  o <- order(prediction)
  deciles <- sort(rep(1:groups,length=length(Y)))

  #OBSERVED AND EXPTECTED
  N <- tapply(Y[o],deciles,length)
  O <- tapply(Y[o],deciles,sum)
  pi <- tapply(prediction[o],deciles,mean)

  num <- (O-N*pi)^2
  denom <- N*pi*(1-pi)

  sum(num/denom)

  chisq = sum(num/denom)

  P = 1 - pchisq(chisq, groups - 2)

  return(list(chisq=chisq,p.value=P))
}

gof.lexpit <- function(object,groups=10){

  Y <- model.frame(object@formula.linear,object@data)[,1]

  #PREDICTED VALUES
  prediction <- predict(object)
  o <- order(prediction)
  deciles <- sort(rep(1:groups,length=length(Y)))

  #OBSERVED AND EXPTECTED
  N <- tapply(Y[o],deciles,length)
  O <- tapply(Y[o],deciles,sum)
  pi <- tapply(prediction[o],deciles,mean)

  num <- (O-N*pi)^2
  denom <- N*pi*(1-pi)

  sum(num/denom)

  chisq = sum(num/denom)

  P = 1 - pchisq(chisq, groups - 2)

  return(list(chisq=chisq,p.value=P))
}


dispersion.blm <- function(object){

     if(class(object)[1]=="lexpit"){
  	X <- model.matrix(object@formula.linear,object@data)
        Z <- model.matrix(object@formula.expit,object@data)
         if(attr(terms(object@formula.linear),"intercept")==1){
           X = X[,-1]
         }
         if(!is.matrix(X)){
           X = matrix(X,ncol=1)
         }
        
        X = cbind(X,Z)
	Y <- model.frame(object@formula.linear,object@data)[,1]
	covariate.class <- apply(X,1,function(x){paste(x,sep="",collapse="")})

        }
     else{
         
	X <- model.matrix(object@formula,object@data)
	Y <- model.frame(object@formula,object@data)[,1]

        }
  
	covariate.class <- apply(X,1,function(x){paste(x,sep="",collapse="")})
	
	pi = predict(object)

	O = tapply(Y,covariate.class,sum)

	N = tapply(Y,covariate.class,length)
	E = N*tapply(pi,covariate.class,mean)
	
	n = length(O)
	p = length(as.vector(coef(object)))
	X2 = sum((O-E)^2/E)
	pearson.df = n-p
	pearson.p = 1-pchisq(X2,df=pearson.df)
	
	D = 2*(sum((Y*log(Y/pi))[Y!=0])+sum(((1-Y)*log((1-Y)/(1-pi)))[Y!=1]))
	deviance.df = sum(N)-p
	deviance.p = 1-pchisq(D,df=deviance.df)

	list(
		observed = O,
		expected = E,
		deviance=D,
		pearson=X2,
		pearson.df=pearson.df,
		deviance.df=deviance.df,
		pearson.p=pearson.p,
		deviance.p=deviance.p
		)

}

setGeneric("gof",function(object){standardGeneric("gof")})

setMethod("gof","blm",dispersion.blm)
setMethod("gof","lexpit",dispersion.blm)


logistic.dispersion <- function(f,data,fit){

	X <- model.matrix(f,data)
        if(!is.matrix(X)) X = matrix(X,ncol=1)
    
	Y <- model.frame(f,data)[,1]
        
	covariate.class <- apply(X,1,function(x){paste(x,sep="",collapse="")})
	
	pi = expit(X%*%coef(fit))

	O = tapply(Y,covariate.class,sum)

	N = tapply(Y,covariate.class,length)
	E = N*tapply(pi,covariate.class,mean)
	
	n = length(O)
	p = length(as.vector(coef(fit)))
	X2 = sum((O-E)^2/E)
	pearson.df = n-p
	pearson.p = 1-pchisq(X2,df=pearson.df)
	
	D = 2*(sum((Y*log(Y/pi))[Y!=0])+sum(((1-Y)*log((1-Y)/(1-pi)))[Y!=1]))
	deviance.df = sum(N)-p
	deviance.p = 1-pchisq(D,df=deviance.df)

	list(
		observed = O,
		expected = E,
		deviance=D,
		pearson=X2,
		pearson.df=pearson.df,
		deviance.df=deviance.df,
		pearson.p=pearson.p,
		deviance.p=deviance.p
		)

}

setMethod("coef","blm",function(object,...) {object@fit$par})
setMethod("coef","lexpit",function(object,...) {object@fit$par})


EO.method <- function(object,alpha=.05,group=NULL,...){

  z = qnorm(1-alpha/2)

  if(class(object)[1]=="blm"){
    f <- object@formula
  }
  else{
    f <- object@formula.linear
  }
  
  Y <- model.frame(f,object@data)[,1]
  pi <- predict(object)

  if(is.null(group)){
   E <- sum(pi)
   O <- sum(Y)
  }
  else{
   E <- tapply(pi,group,sum)
   O <- tapply(Y,group,sum)
  }
  
  lower = E/O*exp(-z*sqrt(1/O))
  upper = E/O*exp(+z*sqrt(1/O))
      
  list(EoverO=E/O,
       lower=lower,
       upper=upper,
       conf.int=(1-alpha)*100
       )
}

setGeneric("EO",function(object,alpha=.05,group=NULL,...){standardGeneric("EO")})
setMethod("EO","blm",EO.method)
setMethod("EO","lexpit",EO.method)

match.constrained <- function(object){ #LIST WITH INDEX OF ACTIVE CONSTRAINTS
  
f.match.constrained <- function(object,active){

  if(is.null(active)) warning("There are no active constraints.")
  
  if(class(object)=="blm"){
    X <- model.matrix(object@formula,object@data)
  }
  else{
    X <- model.matrix(object@formula.linear,object@data)
    if(all(X[,1]==1)){
       X = X[,-1] #REMOVE INTERCEPT TERM
     }
    if(!is.matrix(X)) X = matrix(X,ncol=1)     
    Z <- model.matrix(object@formula.expit,object@data)
    X <- cbind(X,Z)
  }

match.per.index <- function(active){
   row.index <- apply(X,1,function(x){all(x==active)})
   which(row.index)
 }

  if(!is.matrix(active)){
    as.numeric(match.per.index(active))
  }
  else{
   apply(active,1,match.per.index)
  }
}

if(is.matrix(object@active.constraints$active)){
  active.list <- split(object@active.constraints$active,1:nrow(object@active.constraints$active))

  summary.active <- mapply(f.match.constrained,active=active.list,
         MoreArgs=list(object=object),SIMPLIFY=FALSE)

  summary.active$influence <- object@active.constraints$lambda

  summary.active
}
else if(is.constrained(object)){
  active.list <- split(object@active.constraints$active,1)

  summary.active <- mapply(f.match.constrained,active=active.list,
         MoreArgs=list(object=object),SIMPLIFY=FALSE)

  summary.active$influence <- object@active.constraints$lambda

  summary.active
}
else{
  warning("There were no active constraints for the fitted model.")
  NULL
}

}


setGeneric("which.constrained",function(object,active){standardGeneric("which.constrained")})
setMethod("which.constrained","blm",match.constrained)
setMethod("which.constrained","lexpit",match.constrained)

check.constrained <- function(object){
  !any(is.null(object@active.constraints$active))
}

setGeneric("is.constrained",function(object){standardGeneric("is.constrained")})
setMethod("is.constrained","blm",check.constrained)
setMethod("is.constrained","lexpit",check.constrained)

influence.blm <- function(model){

  p <- predict(model)
  X <- model.matrix(model@formula,model@data)
  Y <- model.frame(model@formula,model@data)[,1]
  res <- Y-p  
  nc <- ncol(X)
  a <- 1/(p*(1-p)) 

  b <- b.factor(model)                
  w <- X*(model@weights*as.numeric(a)*as.numeric(res)+b)

  w.bar <- apply(w,1,crossprod)
  w.bar <- w.bar/sum(w.bar)
  names(w.bar) <- 1:nrow(X)
    
  list(influence=w,summary.influence=w.bar)
}


influence.lexpit <- function(model){

  p <- predict(model)
  X <- model.matrix(model@formula.linear,model@data)
  if(all(X[,1]==1)) X = X[,-1] #REMOVE INTERCEPT TERM
  if(!is.matrix(X)) X = matrix(X,ncol=1)     
  Y <- model.frame(model@formula.linear,model@data)[,1]
  Z <- model.matrix(model@formula.expit,model@data)

  res <- Y-p
  px <- ncol(X)
  pz <- ncol(Z)
  nc <- px+pz
  a <- 1/(p*(1-p))
  gamma = coef(model)[(px+1):(px+pz)]
  c <- exp(Z%*%gamma)/(1+exp(Z%*%gamma))^2
  Z <- Z*as.numeric(c)

  XZ <- cbind(X,Z)

  b <- b.factor(model)
  w <- XZ*(model@weights*as.numeric(a)*as.numeric(res)+b)

  w.bar <- apply(w,1,crossprod)
  w.bar <- w.bar/sum(w.bar)
  names(w.bar) <- 1:nrow(X)
   
  list(influence=w,summary.influence=w.bar)
}

setMethod("influence","blm",influence.blm)
setMethod("influence","lexpit",influence.lexpit)
