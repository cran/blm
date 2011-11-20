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
