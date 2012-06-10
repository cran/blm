LRT.blm <- function(object,...){
	
blm.lrt <- function(formula, object,...){
	
	# GIVEN NEW FORMULA AND PROCESSED DATA RETURN LOG-LIKELIHOOD
	
LL <- function(Y,p,w){
		l <- w*(Y*logit(p)+log(1-p))
	-sum(l)
}
	
	Y <- model.frame(formula, object@data)[,1]
	X <- model.matrix(formula, object@data)
	
	# STANDARDIZE WEIGHTS FOR STABILITY
	w <- cbind(object@weights/mean(object@weights))				
	beta.init <- object@par.init[match(colnames(X), names(object@coef))]
	
	fit <- blm.optim(Y,X,w,beta.init,...)
	
-LL(Y,X%*%fit$par,cbind(object@weights))	
}

	warn.option <- getOption("warn")
	options("warn"=-1)

	f <- object@formula	
	labels <- attr(terms(f),"term.labels")
	formulas <- lapply(labels, function(term) update(f, paste("~.-",term, sep="")))
	
	LL <- sapply(formulas, function(new.formula) blm.lrt(new.formula, object,...))
	
	LRTs <- 2*(object@loglik-LL)
	
	table <- cbind(LRT=LRTs, pvalue=1-pchisq(LRTs,df=1))
	row.names(table) <- labels
	
	options("warn"=warn.option)
	
table
}