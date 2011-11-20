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

setGeneric("gof",function(object,groups=10){standardGeneric("gof")})

setMethod("gof","blm",gof.blm)
setMethod("gof","lexpit",gof.lexpit)
