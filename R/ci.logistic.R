ci.logistic.transform <- function(object,C1,C0=NULL,g1=expit,g0=expit,alpha=0.05){

ci.logistic.blm <- function(object,C,alpha=0.05){

  z = qnorm(1 - alpha/2)
   
  fit <- glm(object@formula,object@data,family="binomial")

  est <- coef(fit)
  V <- vcov(fit)
  
  lower.predictor = est-z*sqrt(diag(V))
  upper.predictor = est+z*sqrt(diag(V))

  lower = C%*%lower.predictor
  upper = C%*%upper.predictor

  list(
       lower=lower,
       upper=upper
       )

}

ci.logistic.lexpit <- function(object,C.linear,C.expit,alpha=0.05){

  z = qnorm(1 - alpha/2)

  f.expit <- paste("~.+",paste(object@formula.expit)[3],collapse="")
  f <- update(object@formula.linear,f.expit)

  fit <- glm(f,object@data,family="binomial")

  est <- coef(fit)
  V <- vcov(fit)
  
  lower.predictor = est-z*sqrt(diag(V))
  upper.predictor = est+z*sqrt(diag(V))

  C = c(C.expit[1],C.linear,C.expit[-1])  
  
  lower = C%*%lower.predictor
  upper = C%*%upper.predictor

  list(
       lower=lower,
       upper=upper
       )

}

if(all(is.null(C0))){

  if(class(object)[1]=="blm"){
   C1 <- ci.logistic.blm(object,C1,alpha)
  }
  else{
    C1 <- ci.logistic.lexpit(object,C1$linear,C1$expit,alpha)
   }
  
  lower = g1(C1$lower)
  upper = g1(C1$upper)
   }
else{

  if(class(object)[1]=="blm"){
   C1 <- ci.logistic.blm(object,C1,alpha)
   C0 <- ci.logistic.blm(object,C0,alpha)
  }
  else{
    C1 <- ci.logistic.lexpit(object,C1$linear,C1$expit,alpha)
    C0 <- ci.logistic.lexpit(object,C0$linear,C0$expit,alpha)
   }
  
  lower = g1(C1$lower)-g0(C0$lower)
  upper = g1(C1$upper)-g0(C0$upper)
}

  list(
       lower = lower,
       upper = upper
  )
}

