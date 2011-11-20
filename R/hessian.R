b.factor <- function(object){

  if(all(is.null(object@active.constraints$active))){
    B <- rep(0,nrow(object@data))
  }
  else{
   active.cases <- which.constrained(object)
   active.cases$influence <- NULL
   m <- sapply(active.cases,length) 
   lambdas <- rep(object@active.constraints$lambda,m)
   index <- unlist(active.cases)
   B <- rep(0,nrow(object@data))
   B[index] <- lambdas*rep(1/m,m)

   if(class(object)=="blm"){
     Y <- model.frame(object@formula,object@data)[,1]
   }
   else{
     Y <- model.frame(object@formula.linear,object@data)[,1]
   }
  
   B*ifelse(Y==1,-1,1)           
 }
}


weighted.vcov.blm <- function(object){

  p <- predict(object)
  X <- model.matrix(object@formula,object@data)
  Y <- model.frame(object@formula,object@data)[,1]
  res <- Y-p  
  nc <- ncol(X)
  a <- 1/(p*(1-p)) #WHAT IF INFINITE

  b <- b.factor(object)                
  w <- X*(object@weights*as.numeric(a)*as.numeric(res)+b)

  strata <- 1:nrow(X)
  w.split <- split(as.data.frame(w),strata)
  w.bar <- apply(w,2,mean)

  W2 <- sapply(w.split,function(x){
    y <- apply(x,2,sum)
    outer(y-w.bar,y-w.bar)
  })

  
  matrix(apply(W2,1,sum),nc,nc)   
}


weighted.vcov.lexpit <- function(object){

  p <- predict(object)

  X <- model.matrix(object@formula.linear,object@data)
  if(all(X[,1]==1)) X = X[,-1] #REMOVE INTERCEPT TERM
  if(!is.matrix(X)) X = matrix(X,ncol=1)     
  Y <- model.frame(object@formula.linear,object@data)[,1]
  Z <- model.matrix(object@formula.expit,object@data)

  res <- Y-p
  px <- ncol(X)
  pz <- ncol(Z)
  nc <- px+pz
  a <- 1/(p*(1-p)) #WHAT IF INFINITE?
  gamma = coef(object)[(px+1):(px+pz)]
  c <- exp(Z%*%gamma)/(1+exp(Z%*%gamma))^2
  Z <- Z*as.numeric(c)

  XZ <- cbind(X,Z)
  
  b <- b.factor(object)
  w <- XZ*(object@weights*as.numeric(a)*as.numeric(res)+b)

  strata <- 1:nrow(X)  
  w.split <- split(as.data.frame(w),strata)
  w.bar <- apply(w,2,mean)

  W2 <- sapply(w.split,function(x){
    y <- apply(x,2,sum)
    outer(y-w.bar,y-w.bar)
  })

  
  matrix(apply(W2,1,sum),nc,nc)   
}

hessian.blm <- function(object){

  beta = object@fit$par
  X = model.matrix(object@formula,object@data)
  Y = model.frame(object@formula,object@data)[,1]
  
	x.b = X%*%beta

   h  <- function(i,j){
	 				 sum(-X[,i]*X[,j]*(Y/(x.b)^2+(1-Y)/(1-x.b)^2))
		}
				
	H <- diag(length(beta))

	for(i in 1:length(beta)){
	      for(j in 1:length(beta)){
	      	    if(j>=i) H[i,j] <- h(i,j)
	      }
	}
	
	H[H==0] <- t(H)[H==0]

	H 
}


hessian.lexpit <- function(object){

        n = length(object@fit$par)
  
        X <- model.matrix(object@formula.linear,object@data)
        Y <- model.frame(object@formula.linear,object@data)[,1]
        Z <- model.matrix(object@formula.expit,object@data)
  			
  			if(attr(terms(object@formula.linear),"int")) X = X[,-1]
        p = ifelse(is.matrix(X),ncol(X),1)
        X = matrix(X,ncol=p)
        q = ncol(Z)
        
        dot.expit = function(x){exp(x)/((1+exp(x))^2)}
        ddot.expit = function(x){dot.expit(x)*(1-exp(x))/(1+exp(x))}
        
        beta = object@fit$par[1:p]
        gamma = object@fit$par[(p+1):(p+q)]
        
	      pi = X%*%beta+expit(Z%*%gamma)
             
	      dexpit = dot.expit(Z%*%gamma)
             
              ddexpit = ddot.expit(Z%*%gamma)
             
				
				h.beta <- function(i,j){
	 				 sum(-X[,i]*X[,j]*(Y/(pi)^2+(1-Y)/(1-pi)^2))
				}

      	h.beta.gamma <-  function(i,j){
	       sum(-X[,i]*Z[,j]*dexpit*(Y/(pi)^2+(1-Y)/(1-pi)^2))
	    }
	
	      h.gamma <-  function(i,j){
	 sum(Z[,i]*Z[,j]*ddexpit*(Y/(pi)-(1-Y)/(1-pi))-Z[,i]*Z[,j]*dexpit^2*(Y/(pi)^2+(1-Y)/(1-pi)^2))
	  }
	
	
	H <- diag(n)

	for(i in 1:p){
				
	      for(j in i:n){
	      	    if(j<=p){
	      	      H[i,j] <- h.beta(i,j)
	      	          }
	      	     else{
	      	      H[i,j] <- h.beta.gamma(i,j-p)
	      	      }
	      }
	}
	
	for(i in (p+1):(p+q)){
	
	      for(j in i:(p+q)){
	      	    H[i,j] <- h.gamma(i-p,j-p)
	      }
	}
        	
	H[H==0] <- t(H)[H==0]
        
	H 
}
