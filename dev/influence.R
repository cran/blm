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

vcov.blm.revised <- function(beta,Y,X,w){
	
	x.params <- ncol(X)
	p <- X%*%beta
	A <- 1/nu(p)*(Y-p)
	nu.A <- (A+1)/(nu(p))
	W <- diag(as.numeric(w*nu.A))
	HX <- matrix(-t(X)%*%W%*%X,x.params,x.params)
	
solve(-HX)
}



vcov.lexpit.revised <- function(beta,gamma,Y,X,Z,w){
	
	x.params <- ncol(X)
	z.params <- ncol(Z)
	n.params <- x.params+z.params
	
	p.beta <- X%*%beta
	p.gamma <- expit(Z%*%gamma)
	p <- p.beta+p.gamma
	A <- 1/nu(p)*(Y-p)
	nu.A <- (A+1)/(nu(p))
	
	W.beta <- diag(as.numeric(w*nu.A))
	W.gamma <- diag(as.numeric(w*ddexpit(p.gamma)*A-dexpit(p.gamma)^2*nu.A))
	W.beta.gamma <- diag(as.numeric(w*dexpit(p.gamma)*nu.A))
	
	HX <- -t(X)%*%W.beta%*%X
	HXZ <- -t(X)%*%W.beta.gamma%*%Z
	HZ <- t(Z)%*%W.gamma%*%Z
	
	H <- matrix(0,n.params,n.params)
	
	H[1:x.params,1:x.params] <- HX
	H[1:x.params,(x.params+1):(x.params+z.params)] <- HXZ
	H[(x.params+1):(x.params+z.params),1:x.params] <- t(HXZ)
	H[(x.params+1):(x.params+z.params),(x.params+1):(x.params+z.params)] <- HZ
	
solve(-H)
}

vcov.blm.revised.strata <- function(beta,Y,X,w,strata){

	p <- ncol(X)
	n <- nrow(X)
	H <- matrix(0,p,p)
	
	for(i in levels(strata)){
		n.strata <- sum(strata==i)
		multiplier <- (n-1)/(n-p)*(n.strata/(n.strata-1))
		H <- H+multiplier*solve(vcov.blm.revised(beta,cbind(Y[strata==i,]),cbind(X[strata==i,]),cbind(w[strata==i,])))
	}

solve(H)
}

vcov.lexpit.revised.strata <- function(beta,gamma,Y,X,Z,w,strata){

	p <- ncol(X)+ncol(Z)
	n <- nrow(X)
	H <- matrix(0,p,p)
	
	for(i in levels(strata)){
		n.strata <- sum(strata==i)
		multiplier <- (n-1)/(n-p)*(n.strata/(n.strata-1))
		H <- H+multiplier*solve(vcov.lexpit.revised(beta,gamma,
					cbind(Y[strata==i,]),cbind(X[strata==i,]),cbind(Z[strata==i,]),cbind(w[strata==i,])))
	}

solve(H)
}

deviates.lexpit.revised <- function(beta,gamma,Y,X,Z,w){
	

	p.beta <- X%*%beta
	p.gamma <- expit(Z%*%gamma)
	p <- p.beta+p.gamma
	
	# WEIGHTED SECOND DERIVATIVES
	XX <- apply(X,1,function(x)outer(x,x))
	if(!is.matrix(XX)) XX <- rbind(XX)
	XX <- XX*matrix(w,nr=nrow(XX),ncol=ncol(XX),byrow=T)
	HX <- matrix(rowSums(XX),ncol(X),ncol(X))
	
	ZZ <- apply(Z,1,function(x)outer(x,x))
	if(!is.matrix(ZZ)) ZZ <- rbind(ZZ)
	ZZ <- ZZ*matrix(w*p.gamma*(1-p.gamma),nr=nrow(ZZ),ncol=ncol(ZZ),byrow=T)
	HZ <- matrix(rowSums(ZZ),ncol(Z),ncol(Z))
	
	x.score <- w*(Y-p.beta)
	z.score <- w*(Y-p.gamma)
	
	X.score <- X*matrix(x.score,nr=nrow(X),nc=ncol(X))
	Z.score <- Z*matrix(z.score,nr=nrow(Z),nc=ncol(Z))
	
	X.score <- apply(X.score,1,function(x) solve(HX)%*%x)
	Z.score <- apply(Z.score,1,function(x) solve(HZ)%*%x)
	
	if(!is.matrix(X.score)) 
		X.score <- cbind(X.score)
	else
		X.score <- t(X.score)
	
	if(!is.matrix(Z.score)) 
		Z.score <- cbind(Z.score)
	else
		Z.score <- t(Z.score)
			
	list(
		linear = X.score,
		expit = Z.score
		)
}



influence.lexpit.revised <- function(beta.deviates,gamma.deviates,beta,gamma,strata){

	V <- cbind(beta.deviates,gamma.deviates)	
	V.split <- split(data.frame(V), strata)
	
	f <- function(x){
		n <- nrow(x)
		means <- colSums(x)/n
		X <- apply(x,1,function(y)outer(y-means,y-means))
		X <- X*n/(n-1)
		if(is.matrix(X))
			rowSums(X)
		else
			sum(X)
	}
			
	V <- sapply(V.split,f)
		
	p <- ncol(beta.deviates)
	q <- ncol(gamma.deviates)
	
	V <- matrix(rowSums(V),p+q,p+q)	
	vcov.beta <- matrix(V[1:p,1:p],p,p)
	vcov.gamma <- matrix(V[(p+1):(p+q),(p+1):(p+q)],q,q)

list(
	beta = vcov.beta,
	gamma = vcov.gamma
)
}

influence <- function(beta.deviates,gamma.deviates,strata){
	
	# GROUP BY STRATA
	split.beta <- split(data.frame(t(beta.deviates)),strata)
	split.gamma <- split(data.frame(t(gamma.deviates)),strata)
	
	p <- nrow(beta.deviates)
	q <- nrow(gamma.deviates)
	
	f <- function(x){
		n <- nrow(x)
		means <- colSums(x)/n
		X <- apply(x,1,function(y)outer(y-means,y-means))
		X <- X*n/(n-1)
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


influence.old <- function(beta.deviates,gamma.deviates,strata){
	
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
		
	p <- nrow(beta.deviates) # NUMBER OF PARAMETERS
	
	f <- function(x){ 
			# RETURNS SQUARED INFLUENCE BY STRATA
			# EACH STRATUM IS A COLUMN
		n <- nrow(x)
		means <- colSums(x)/n
		X <- apply(x,1,function(y)outer(y-means,y-means))
		X <- X*n/(n-1)
		if(is.matrix(X))
			rowSums(X)
		else
			sum(X)
	}
	
		
	V.beta <- sapply(split.beta,f)
	if(!is.matrix(V.beta)) V.beta <- rbind(V.beta)
	
matrix(rowSums(V.beta),p,p)
}



blm.old.influence <- function(beta.deviates,strata){
	
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


