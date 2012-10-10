influence.blm <- function(formula, data, weights=NULL){
	fit <- do.call("glm", args=list(formula=formula, data=data, 
												weights=weights,family=binomial(link=make.link("identity"))))
	lm.influence(fit)$coefficients
}

influence.lexpit <- function(formula, data, weights=NULL){
	
	fit <- do.call("glm", args=list(formula=formula, data=data, 
												weights=weights,family="binomial"))
	lm.influence(fit)$coefficients
}

vcov.influence.blm <- function(formula, data){
	influence <- influence.blm(formula, data)
t(influence)%*%influence
}

vcov.influence.lexpit <- function(formula.linear, formula.expit, data){

	influence.linear <- influence.blm(formula.linear, data)
	influence.expit <- influence.lexpit(formula.expit, data)
	
	influence <- cbind(influence.linear[,-1],influence.expit)
	
t(influence)%*%influence
}


vcov.influence.blm.strata <- function(formula, data, weights, strata){
	
	influence <- influence.blm(formula, data, weights)
	means <- influence
	size <- table(strata)[strata]
	size <- size/(size-1)
	
	for(i in 1:ncol(influence)){
		means[,i] <- (tapply(influence[,i], strata, mean)[strata])
		influence[,i] <- (influence[,i]-means[,i])*size
	}
	
t(influence)%*%influence
}


vcov.influence.lexpit.strata <- function(formula.linear, formula.expit, data, weights, strata){

	influence.linear <- influence.blm(formula.linear, data, weights)
	influence.expit <- influence.lexpit(formula.expit, data, weights)
	influence <- cbind(influence.linear[,-1],influence.expit)
	
	means <- influence
	size <- table(strata)[strata]
	size <- size/(size-1)
	
	for(i in 1:ncol(influence)){
		means[,i] <- (tapply(influence[,i], strata, mean)[strata])
		influence[,i] <- (influence[,i]-means[,i])*size
	}
	
t(influence)%*%influence
}

vcov.blm.big <- function(formula, data, weights=NULL){
	fit <- do.call("glm", args=list(formula=formula, data=data, 
												weights=weights,family=binomial(link=make.link("identity"))))
vcov(fit)
}

vcov.lexpit.big <- function(formula.linear, formula.expit, data, weights=NULL){

	fit <- do.call("glm", args=list(formula=formula.linear, data=data, 
												weights=weights,family=binomial(link=make.link("identity"))))

   vcov.linear <- vcov(fit)
   
   	fit <- do.call("glm", args=list(formula=formula.expit, data=data, 
												weights=weights,family="binomial"))
	vcov.expit <- vcov(fit)

	p <- nrow(vcov.linear)	
	q <- nrow(vcov.expit)
	
	V <- matrix(0, p+q-1, p+q-1)
	
	V[(1:(p-1)),(1:(p-1))] <- vcov.linear[2:p,2:p] 
	V[p:(p+q-1),p:(p+q-1)] <- vcov.expit
	
V
}


