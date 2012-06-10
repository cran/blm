predict.blm <- function(object){
	
	if(class(object)=="lexpit"){
		X <- model.matrix(object@formula.linear,object@data)
		Z <- model.matrix(object@formula.expit,object@data)
	
		x.has.intercept <- attr(terms(object@formula.linear),"intercept")==1
	
		if(x.has.intercept) X <- X[,-1]
	
		if(!is.matrix(X)) X <- cbind(X)
		if(!is.matrix(Z)) Z <- cbind(Z)

		p <- X%*%object@coef.linear+expit(Z%*%object@coef.expit)
	}
	else{
		X <- model.matrix(object@formula,object@data)
		p <- X%*%object@coef
	}
p
}

setMethod("predict","blm",predict.blm)
setMethod("predict","lexpit",predict.blm)