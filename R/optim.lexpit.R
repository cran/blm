# GIVEN GAMMA COMPUTE BETA, GAMMA, LOGLIK

optim.lexpit <- function(beta.init,gamma.init,Y,X,Z,w,...){
	
	# NEED TO MATCH ARGUMENTS FOR OPTIM/NLM
	
	LL <- function(Y,p,w){
		l <- w*(Y*logit(p)+log(1-p))
	-sum(l)
}

	beta.offset <- expit(Z%*%gamma.init)
	beta.optim <- stage1.optim(Y,X,w,beta.offset,beta.init)
		
	gamma.offset <- X%*%beta.optim$par
	gamma.optim <- stage2.optim(Y,Z,gamma.offset,w,gamma.init)$est
		
	l <- LL(Y,gamma.offset+expit(Z%*%gamma.optim),w)
	
list(beta=beta.optim$par,gamma=gamma.optim,loglik=l,barrier.value=beta.optim$barrier.value)	
}
