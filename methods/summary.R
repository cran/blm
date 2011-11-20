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
