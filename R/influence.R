deviates <- function(beta,gamma,Y,X,Z,w){
	
	p.beta <- X%*%beta
	p.gamma <- expit(Z%*%gamma)
	p <- p.beta+p.gamma
	
	C1 <- w*dexpit(Z%*%gamma)^2*(Y*1/(p*(1-p))^2*(1-2*p)+1/(1-p)^2)
	C2 <- w*ddexpit(Z%*%gamma)*(Y*1/(p*(1-p))-1/(1-p))
	ZZ <- apply(Z,1,function(z)outer(z,z))
	C <- matrix(C2-C1,ncol(Z)^2,nrow(Z),byrow=T)
	
	H.gamma <- matrix(rowSums(C*ZZ),ncol(Z),ncol(Z))
	H.gamma <- solve(-H.gamma)
	
	C <- -w*(Y*1/(p*(1-p))^2*(1-2*p)+1/(1-p)^2)
	XX <- apply(X,1,function(z)outer(z,z))
	C <- matrix(C,ncol(X)^2,nrow(X),byrow=T)
	
	H.beta <- matrix(rowSums(C*XX),ncol(X),ncol(X))
	H.beta <- solve(-H.beta)
	
	beta.score <- w*(Y*1/(p*(1-p))-1/(1-p))
	gamma.score <- beta.score*dexpit(Z%*%gamma)

	HX <- apply(X,1,function(x)H.beta%*%x)
	HZ <- apply(Z,1,function(x)H.gamma%*%x)

	beta.deviates <- HX*matrix(beta.score,ncol(X),nrow(X),byrow=T)
	gamma.deviates <- HZ*matrix(gamma.score,ncol(Z),nrow(Z),byrow=T)
	
	list(
		linear = beta.deviates,
		expit = gamma.deviates
	)
}

influence <- function(beta.deviates,gamma.deviates,strata){
	
	# GROUP BY STRATA
	split.beta <- split(data.frame(t(beta.deviates)),strata)
	split.gamma <- split(data.frame(t(gamma.deviates)),strata)
	
	p <- nrow(beta.deviates)
	q <- nrow(gamma.deviates)
	
	f <- function(x){
		means <- colSums(x)/nrow(x)
		X <- apply(x,1,function(y)outer(y-means,y-means))
		if(is.matrix(X))
			rowSums(X)
		else
			sum(X)
	}
	
		
	V.beta <- sapply(split.beta,f)
	V.gamma <- sapply(split.gamma,f)
	
	if(!is.matrix(V.beta)) V.beta <- rbind(V.beta)
	if(!is.matrix(V.gamma)) V.gamma <- rbind(V.gamma)
	
	n <- length(split.beta)
	
	vcov.beta <- matrix(rowSums(V.beta),p,p)
	vcov.gamma <- matrix(rowSums(V.gamma),q,q)


list(
	beta = vcov.beta,
	gamma = vcov.gamma
)
}


blm.deviates <- function(beta,Y,X,w){
	
	p <- X%*%beta
	
	C <- -w*(Y*1/(p*(1-p))^2*(1-2*p)+1/(1-p)^2)
	XX <- apply(X,1,function(z)outer(z,z))
	C <- matrix(C,ncol(X)^2,nrow(X),byrow=T)
	
	H.beta <- matrix(rowSums(C*XX),ncol(X),ncol(X))
	H.beta <- solve(-H.beta)
	
	beta.score <- w*(Y*1/(p*(1-p))-1/(1-p))
	HX <- apply(X,1,function(x)H.beta%*%x)

HX*matrix(beta.score,ncol(X),nrow(X),byrow=T)
}

blm.influence <- function(beta.deviates,strata){
	
	# GROUP BY STRATA
	split.beta <- split(data.frame(t(beta.deviates)),strata)
	p <- nrow(beta.deviates)
	
	f <- function(x){
		means <- colSums(x)/nrow(x)
		X <- apply(x,1,function(y)outer(y-means,y-means))
		if(is.matrix(X))
			rowSums(X)
		else
			sum(X)
	}
	
		
	V.beta <- sapply(split.beta,f)

	if(!is.matrix(V.beta)) V.beta <- rbind(V.beta)
	
	n <- length(split.beta)
	
matrix(rowSums(V.beta),p,p)
}


