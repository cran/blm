starting.values.lexpit <- function(f.linear,f.expit,data,par.init)
{

        
     	X <- model.matrix(f.linear,data)
        X.names <- colnames(X)

        if(all(X[,1]==1)){
          X = X[,-1] #REMOVE INTERCEPT TERM
          X.names = X.names[-1]
        }
        if(!is.matrix(X)) X = matrix(X,ncol=1)     
	Z <- model.matrix(f.expit,data)
        Z.names <- colnames(Z)
        
        p = ncol(X)
        q = ncol(Z)
        not.feasible <- FALSE
                  
         if(!missing(par.init)){
           beta = matrix(par.init$linear,ncol=1)
           gamma = matrix(par.init$expit,ncol=1)
           not.feasible = any((X%*%beta+expit(Z%*%gamma)<=0)|(X%*%beta+expit(Z%*%gamma)>=1))
           }

         if(missing(par.init)|not.feasible){
           if(not.feasible){
             warning("Initial values not feasible. Selecting alternative.")
           }

            gamma <- coef(glm(f.expit,data,family="binomial"))
            pre <- expit(Z%*%gamma)
            Y <- model.matrix(f.linear,data)[,1]-pre
            beta <- coef(glm(f.linear,data,family="gaussian"))
            if(length(beta)>p) beta = beta[-1]
            not.feasible = any((X%*%beta+expit(Z%*%gamma)<=0)|(X%*%beta+expit(Z%*%gamma)>=1))
          
           while(not.feasible){
             beta = beta*runif(1,.9,1)
             not.feasible = any((X%*%beta+expit(Z%*%gamma)<=0)|(X%*%beta+expit(Z%*%gamma)>=1))
           }
         }
         
              return(list(par.start=c(beta,gamma),
               names = c(X.names,Z.names),
               not.feasible=not.feasible))

}

starting.values.blm <- function(f,data,par.init)
{

         X <- model.matrix(f,data)
         p <- ncol(X)
         if(!is.matrix(X)) X <- matrix(X,ncol=1)

         A <- unique(X)
         not.feasible <- FALSE
                  
         if(!missing(par.init)){
           beta = matrix(par.init,ncol=1)
           not.feasible = any((A%*%beta<=0)|(A%*%beta>=1))
           }

         if(missing(par.init)|not.feasible){
           if(not.feasible){
             warning("Initial values not feasible. Selecting alternative.")
           }

            beta <- coef(lm(f,data))
       	    not.feasible = any((A%*%beta<=0)|(A%*%beta>=1))
            
            if(not.feasible){
              if(beta[1]<=0|beta[1]>=1) beta[1] = coef(lm(update(f,.~1),data))
            }
           
           while(not.feasible){
             beta[2:p] = beta[2:p]*runif(1,.9,1)
             not.feasible = any((A%*%beta<=0)|(A%*%beta>=1))
           }
         }
         
	 return(list(par.start=beta,
                     names=colnames(A),
                     not.feasible=not.feasible
                    ))
}
