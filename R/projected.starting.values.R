starting.values.blm <- function(f,data,par.init)
{

	X <- model.matrix(f,data)
        X.names <- colnames(X)
        if(!is.matrix(X)) X = matrix(X,ncol=1)
        
        A <- blm.ineq.matrix(f,data)
  
	Y <- model.frame(f,data)[,1]

        beta <- solve(t(X)%*%X)%*%t(X)%*%Y

         if(!missing(par.init)){
           beta = matrix(par.init,ncol=1)
           rownames(beta) = X.names
         }

        A = unique(A)
        
	not.feasible = any((A%*%beta<=0)|(A%*%beta>=1))

	 if(not.feasible){
	   if(!missing(par.init)){
             warning("Initial values were not within feasible region. Feasible starting values selected.")
             beta <- solve(t(X)%*%X)%*%t(X)%*%Y
           }
           
           Amat = rbind(rbind(X,A),-rbind(X,A))
           min <- abs(min(rbind(A,X)%*%beta))
           B = rep(c(min,-1+min),each=(nrow(A)+nrow(X)))
        
           beta <- projectLinear(beta,Amat,B,0)
           not.feasible = any((A%*%beta<=0)|(A%*%beta>=1))
        }

	 return(list(par.start=beta,
                    names=X.names,
                    not.feasible=not.feasible
                    ))
}


starting.values.lexpit <- function(f.linear,f.expit,data,par.init)
{
	
	X <- model.matrix(f.linear,data)
	Y <- model.frame(f.linear,data)[,1]
	Z <- model.matrix(f.expit,data)
        
        X.names <- colnames(X)
        Z.names <- colnames(Z)
        
        int = attr(terms(f.linear),"intercept")
        
	if(!int){
            if(!is.matrix(X)) X = matrix(X,ncol=1)
            X = cbind(rep(1,length(Y)),X)
          }
        else{
          X.names = X.names[-1]
        }
        
        p = ncol(X)-1
        q = ncol(Z)
        

        beta <- solve(t(X)%*%X)%*%t(X)%*%Y
        gamma <- c(log(mean(Y)/(1-mean(Y))),rep(0,q-1))
      	pi.gamma <- expit(gamma[1])
	    
        beta = beta[-1]
	X <- X[,-1]
	if(p==1) X = matrix(X,ncol=1)
				
        if(!missing(par.init)){
            beta <- par.init$linear
            gamma <- par.init$expit
            names(beta) <- X.names
            names(gamma) <- Z.names
         }
	
	A = lexpit.ineq.matrix(f.linear,data)
	A = rbind(X,A)
	A = unique(A)

	not.feasible = any((A%*%beta<=-pi.gamma)|(A%*%beta>=1-pi.gamma))

	if(not.feasible){
	   if(!missing(par.init)){
             warning("Initial values were not within feasible region. Feasible starting values selected.")
             beta <- rep(0,length(beta))
           }
           
           Amat = rbind(A,-A)
           min = 1/1000
           B = rep(c(min-pi.gamma,-1+min+pi.gamma),each=nrow(A))
        
           beta <- projectLinear(beta,Amat,B,0)
           not.feasible = any((A%*%beta<=-pi.gamma)|(A%*%beta>=1-pi.gamma))
 		}
  				    
   return(list(par.start=c(beta,gamma),
               names = c(X.names,Z.names),
               not.feasible=not.feasible))
}
