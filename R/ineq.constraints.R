blm.ineq.matrix <- function(f,data){
  
   X <- model.matrix(f,data)    
   unique(X)
}


blm.constraints <- function(f,data,ineq.mat=NULL){

  # JACOBIAN IS THE FIRST DERIVATIVE
  
  if(is.null(ineq.mat)){
    ineq.mat <- blm.ineq.matrix(f,data)
  }
  
  ineq <- function(par){c(ineq.mat%*%par,1-ineq.mat%*%par)}

  ineq.jac <- function(par){
    
     rbind(ineq.mat,-ineq.mat)
  }

  return(list(ineq.mat=ineq.mat,ineq=ineq,ineq.jac=ineq.jac))
}


lexpit.ineq.matrix <- function(f,data){
		
        A <- model.matrix(f,data)
        
        if(attr(terms(f),"int")) A = A[,-1]
        if(!is.matrix(A)) A = matrix(A,ncol=1)
   A
}


lexpit.constraints <- function(f.linear,f.expit,data,ineq.mat=NULL){

    A <- lexpit.ineq.matrix(f.linear,data)
    B <- model.matrix(f.expit,data)
    p = ncol(A)
    q = ncol(B)
    
  if(is.null(ineq.mat)){
    
    U <- cbind(A,B) #REDUCE TO UNIQUE COVARIATE CLASSES
    
    ineq.mat = list(
      A = matrix(U[,1:p],ncol=p),
      B = U[,(p+1):(p+q)]
      )
  }
  
  ineq <- function(par){
    
    beta = par[1:p]
    gamma = par[(p+1):(p+q)]
    
    c(ineq.mat$A%*%beta+expit(ineq.mat$B%*%gamma),
      1-ineq.mat$A%*%beta-expit(ineq.mat$B%*%gamma))
  }

  ineq.jac <- function(par){

    dot.expit = function(x) {
        exp(x)/(1 + exp(x))^2
    }
     
     beta = par[1:p]
     gamma = par[(p+1):(p+q)]

     gamma.factor <- dot.expit(ineq.mat$B%*%gamma)
     gamma.factor <- matrix(gamma.factor,nrow=nrow(gamma.factor),ncol=q)
                            
     gamma.jac <- rbind(ineq.mat$B*gamma.factor,-ineq.mat$B*gamma.factor)
     beta.jac = rbind(ineq.mat$A,-ineq.mat$A)
     cbind(beta.jac,gamma.jac)    
  }

  return(list(ineq.mat=ineq.mat,ineq=ineq,ineq.jac=ineq.jac))

}

