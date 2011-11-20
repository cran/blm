efron <- function(obs,pred){
  1-sum((obs-pred)^2)/sum((obs-mean(obs))^2)
}

mcfadden <- function(loglik,loglik.null){
  1-loglik/loglik.null
}

mcfadden.adj <- function(loglik,loglik.null,num.params){
  1-(loglik-num.params)/loglik.null
}

cox.snell <- function(loglik,loglik.null,N){
  1-exp(2/N*(loglik.null-loglik))
}

cox.snell.adj <- function(loglik,loglik.null,N){
  cox.snell(loglik,loglik.null,N)/(1-exp(loglik.null)^(2/N))
}

