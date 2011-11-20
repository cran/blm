# NEGATIVE OF LIKELIHOOD EQUATIONS FOR MINIMIZATION OPTIMIZER

blm.loglik <- function(f,data,w){

        X <- model.matrix(f,data)
        Y <- model.frame(f,data)[,1]

        LL <- function(beta){
          
	x.b = X%*%beta

	ll1 <- sum(w*Y*(log(x.b)-log(1-x.b)))
	ll2 <- sum(w*log(1-x.b))

	-(ll1+ll2)
      }
        
    LL
}


blm.dot.loglik <- function(f,data,w){

        X <- model.matrix(f,data)
        Y <- model.frame(f,data)[,1]

        dot.LL <- function(beta){
	x.b = X%*%beta

	g <- sapply(1:ncol(X),function(x){
	  sum(w*Y*X[,x]/(x.b))-sum(w*(1-Y)*X[,x]/(1-x.b))
	})

	-g
      }

        dot.LL
}

lexpit.loglik <- function(f.linear,f.expit,data,w){

	X <- model.matrix(f.linear,data)
        if(all(X[,1]==1)) X = X[,-1] #REMOVE INTERCEPT TERM
        if(!is.matrix(X)) X = matrix(X,ncol=1)     
        Y <- model.frame(f.linear,data)[,1]
        Z <- model.matrix(f.expit,data)

        p = ncol(X)
        q = ncol(Z)

        
        LL <- function(par){
          
        beta = par[1:p]
        gamma = par[(p+1):(p+q)]

 	pi = X%*%beta+expit(Z%*%gamma)

	ll1 <- sum(w*Y*(log(pi)-log(1-pi)))
	ll2 <- sum(w*log(1-pi))

	-(ll1+ll2)

      }

        LL
}

lexpit.dot.loglik <- function(f.linear,f.expit,data,w){

  	X <- model.matrix(f.linear,data)
        if(all(X[,1]==1)) X = X[,-1] #REMOVE INTERCEPT TERM
        if(!is.matrix(X)) X = matrix(X,ncol=1)     
        Y <- model.frame(f.linear,data)[,1]
        Z <- model.matrix(f.expit,data)

        p = ncol(X)
        q = ncol(Z)
        
        dot.LL <- function(par){
          
        dot.expit = function(x){exp(x)/(1+exp(x))^2}
        
        beta = par[1:p]
        gamma = par[(p+1):(p+q)]
        
	      pi = X%*%beta+expit(Z%*%gamma)

	g.beta <- sapply(1:ncol(X),function(x){
	  sum(w*Y*X[,x]/(pi))-sum(w*(1-Y)*X[,x]/(1-pi))
	})

  	g.gamma <- sapply(1:ncol(Z),function(x){
	  sum(w*Y*Z[,x]*dot.expit(Z%*%gamma)/(pi))-sum(w*(1-Y)*Z[,x]*dot.expit(Z%*%gamma)/(1-pi))
	})
  
        -c(g.beta,g.gamma)
      }
        dot.LL
}

