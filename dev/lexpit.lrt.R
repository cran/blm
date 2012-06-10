LRT.lexpit <- function(object,...){
	
lexpit.lrt <- function(formula.linear, formula.expit, object,...){

LL <- function(Y,p,w){
		l <- w*(Y*logit(p)+log(1-p))
	-sum(l)
}

	if(is.null(formula.linear)){	
		Y <- model.frame(formula.expit,object@data)[,1]
		Z <- model.matrix(formula.expit,object@data)
		X <- NULL
		if(!is.matrix(Z)) Z <- cbind(Z)
		beta.init <- NULL
		gamma.init <- object@par.init[(object@p+1):(object@p+object@q)]
	}
	else if(is.null(formula.expit)){
		Y <- model.frame(formula.linear,object@data)[,1]
		X <- model.matrix(formula.linear,object@data)
		Z <- NULL
		if(!is.matrix(X)) X <- cbind(X)
		beta.init <- object@par.init[1:object@p] # ADD INTERCEPT
		beta.init <- c(sum(Y*object@weights)/sum(object@weights),beta.init)
		gamma.init <- NULL
	}
	else{
		
			Y <- model.frame(formula.linear,object@data)[,1]
			X <- model.matrix(formula.linear,object@data)
			Z <- model.matrix(formula.expit,object@data)
	
			x.has.intercept <- attr(terms(formula.linear),"intercept")==1
			x.labels <- attr(terms(formula.linear),"term.labels")
			z.labels <- c("(Intercept)",attr(terms(formula.expit),"term.labels"))
			
			if(x.has.intercept) X <- X[,-1]
			if(!is.matrix(X)) X <- cbind(X)
			if(!is.matrix(Z)) Z <- cbind(Z)
			beta.init <- object@par.init[1:object@p]
			beta.init <- beta.init[match(x.labels, names(object@coef.linear))]
			
			gamma.init <- object@par.init[(object@p+1):(object@p+object@q)]
			gamma.init <- gamma.init[match(z.labels, names(object@coef.expit))]
	}
	
	
	# STANDARDIZE WEIGHTS FOR STABILITY
	w <- cbind(object@weights/mean(object@weights))
	
	i <- 0
	threshold.met <- FALSE
	criterion <- 0
	
	beta <- beta.init
	gamma <- gamma.init
	
	while(i<object@control.lexpit$max.iter&!threshold.met){
		
		fit <- optim.lrt.lexpit(beta,gamma,Y,X,Z,w)
		beta <- fit$beta
		gamma <- fit$gamma
		threshold.met <- abs(fit$loglik-criterion)<object@control.lexpit$tol
		criterion <- fit$loglik
		i <- i+1
	}
	
	if(is.null(formula.linear)){
		-LL(Y,expit(Z%*%gamma),cbind(object@weights))
	}
	else if(is.null(formula.expit)){
		-LL(Y,X%*%beta,cbind(object@weights))
	}
	else{
		-LL(Y,X%*%beta+expit(Z%*%gamma),cbind(object@weights))
	}
}

	warn.option <- getOption("warn")
	options("warn"=-1)
	
	f.linear <- object@formula.linear	
	labels.linear <- attr(terms(f.linear),"term.labels")
	
	f.expit <- object@formula.expit
	labels.expit <- attr(terms(f.expit),"term.labels")
	
	if(length(labels.linear)==1&length(labels.expit)==1){
		formulas.linear <- list(NULL)
		formulas.expit <- list(NULL)
	}
	else if(length(labels.linear)==1){
		formulas.linear <- list(NULL)
		formulas.expit <- lapply(labels.expit, function(term) update(f.expit, paste("~.-",term, sep="")))
	}
	else if(length(labels.expit)==1){
		formulas.expit <- list(NULL)
		formulas.linear <- lapply(labels.linear, function(term) update(f.linear, paste("~.-",term, sep="")))
	}
	else{	
		formulas.expit <- lapply(labels.expit, function(term) update(f.expit, paste("~.-",term, sep="")))
		formulas.linear <- lapply(labels.linear, function(term) update(f.linear, paste("~.-",term, sep="")))
	}
	
	LL.linear <- sapply(formulas.linear, function(new.formula) lexpit.lrt(new.formula, object@formula.expit, object,...))
	LL.expit <- sapply(formulas.expit, function(new.formula) lexpit.lrt(object@formula.linear, new.formula, object,...))
	
	LRTs <- 2*(object@loglik-c(LL.linear, LL.expit))	
	table <- cbind(LRT=LRTs, pvalue=1-pchisq(LRTs,df=1))
	row.names(table) <- c(labels.linear, labels.expit)
	options("warn"=warn.option)
	
table
}