blm <- function(formula,data,na.action=na.omit,
							weights=NULL,strata=NULL,par.init=NULL,...){

na.lexpit <- function(f,data,FUN){

	keep.data <- subset(data,select=all.vars(f))
	keep.data <- FUN(keep.data)
	kept <- match(row.names(keep.data),row.names(data))

list(data = keep.data, kept = kept)
}

LL <- function(Y,p,w){
		l <- w*(Y*logit(p)+log(1-p))
	-sum(l)
}

	data <- na.lexpit(formula,data,FUN=na.action)
	which.kept <- data$kept
	data <- data$data
	
	Y <- model.frame(formula,data)[,1]
	X <- model.matrix(formula,data)
	
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
	
	if(is.null(par.init)){
		beta.init <- rep(0,ncol(X))
		beta.init[1] <- sum(Y*w)/sum(w)
	}
	else{
		beta.init <- par.init
		}
		
	fit <- blm.optim(Y,X,w,beta.init,...)
	beta <- fit$par
		
	devs <- blm.deviates(beta,Y,X,cbind(w))
	vcov <- blm.influence(devs,strata)
	
	devs <- t(devs)
	
	if(attr(terms(formula),"intercept")==1){
		devs <- cbind(devs[,-1])
		colnames(devs) <- colnames(X)[-1]
	}
	else{
		colnames(devs) <- colnames(X)
		}
		
	names(beta) <- colnames(X)
	
	# NULL MODEL, ESTIMATE IS THE OVERALL MEAN
	ll.null <- -LL(Y,sum(Y*w)/sum(w),cbind(weights))
	
	new("blm",
		coef = beta,
		vcov = vcov,
		formula= formula,
		df.residual = nrow(X)-ncol(X),
		data = data,
		which.kept = which.kept,
		y = Y,
		weights = as.numeric(weights),
		strata = strata,
		converged = fit$convergence==0,
		par.init = beta.init,
		loglik = -LL(Y,X%*%beta,cbind(weights)),
		loglik.null = ll.null,
		barrier.value = fit$barrier.value,
		dBeta = devs
	)
	
}

