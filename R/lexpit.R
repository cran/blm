lexpit <- function(formula.linear,formula.expit,data,na.action=na.omit,
					weights=NULL,strata=NULL, par.init = NULL,
					control.lexpit=list(max.iter=1000,tol=1E-7),...){

na.lexpit <- function(f1,f2,data,FUN){

	keep.data <- subset(data,select=unique(c(all.vars(f1),all.vars(f2))))
	keep.data <- FUN(keep.data)
	kept <- match(row.names(keep.data),row.names(data))

list(data = keep.data, kept = kept)
}

LL <- function(Y,p,w){
		l <- w*(Y*logit(p)+log(1-p))
	-sum(l)
}
	data <- na.lexpit(formula.linear,formula.expit,data,FUN=na.action)
	which.kept <- data$kept
	data <- data$data
	
	Y <- model.frame(formula.linear,data)[,1]
	X <- model.matrix(formula.linear,data)
	Z <- model.matrix(formula.expit,data)
	
	x.has.intercept <- attr(terms(formula.linear),"intercept")==1
	
	if(x.has.intercept) X <- X[,-1]
	
	if(!is.matrix(X)) X <- cbind(X)
	if(!is.matrix(Z)) Z <- cbind(Z)
	if(is.null(weights)){
			weights <- rep(1,nrow(X))
			}
	else{
			weights <- weights[which.kept]
			}
	
	# STANDARDIZE WEIGHTS FOR STABILITY
	w <- cbind(weights/mean(weights))
		
	if(is.null(strata)){
			strata <- rep(1,nrow(X))
			}
		else{
			strata <- strata[which.kept]
			}
			
	if(!class(strata)=="factor") strata <- factor(strata)
	
	if(!is.list(par.init)){
		beta.init <- rep(0,ncol(X))
		gamma.init <- do.call(glm,
					list(
						formula=formula.expit,
						data=data,
						family="binomial",
						weight=weights))$coef
						}
	else{
		beta.init <- par.init$linear
		gamma.init <- par.init$expit
	}
	
	i <- 0
	threshold.met <- FALSE
	criterion <- 0
	
	beta <- beta.init
	gamma <- gamma.init
	
	while(i<control.lexpit$max.iter&!threshold.met){
		
		fit <- optim.lexpit(beta,gamma,Y,X,Z,w)
		beta <- fit$beta
		gamma <- fit$gamma
		threshold.met <- abs(fit$loglik-criterion)<control.lexpit$tol
		criterion <- fit$loglik
		i <- i+1
	}
	
	
	coef.deviates <- deviates(beta,gamma,Y,X,Z,cbind(w)) 
	vcov <- influence(coef.deviates[[1]],coef.deviates[[2]],strata)
	
	coef.deviates[[1]] <- t(coef.deviates[[1]])
	colnames(coef.deviates[[1]]) <- colnames(X)
	
	coef.deviates[[2]] <- t(coef.deviates[[2]])
	coef.deviates[[2]] <- cbind(coef.deviates[[2]][,-1])
	colnames(coef.deviates[[2]]) <- colnames(Z)[-1]
	
	names(beta) <- colnames(X)
	names(gamma) <- colnames(Z)
	
	if(ncol(X)==1) names(beta) <- all.vars(formula.linear)[2]
	if(ncol(Z)==1) names(gamma) <- "(Intercept)"
	
	# NULL MODEL, ESTIMATE IS THE OVERALL MEAN
	ll.null <- -LL(Y,sum(Y*w)/sum(w),cbind(weights))
	
	new("lexpit",
		coef.linear = beta,
		coef.expit = gamma,
		vcov.linear = vcov[[1]],
		vcov.expit = vcov[[2]],
		formula.linear = formula.linear,
		formula.expit = formula.expit,
		df.residual = nrow(X)-ncol(X)-ncol(Z),
		p = ncol(X),
		q = ncol(Z),
		data = data,
		which.kept = which.kept,
		y = Y,
		weights = as.numeric(weights),
		strata = strata,
		converged = i<control.lexpit$max.iter,
		par.init = c(beta.init, gamma.init),
		loglik = -LL(Y,X%*%beta+expit(Z%*%gamma),cbind(weights)),
		loglik.null = ll.null,
		barrier.value = fit$barrier.value,
		dBeta = coef.deviates
	)
	
}

