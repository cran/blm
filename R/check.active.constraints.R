check.active.constraints <- function(object,tol=1/1000){

  constraints <- object@constraints$ineq(coef(object))

  if(any(abs(constraints)<tol)){

    i <- which(abs(constraints)<tol)
    A <- object@constraints$ineq.jac(coef(object))
    
    list(active.constraints=object@ineq[i,],boundary.value=constraints[i])
  }
  else{
    list(active.constraints=NULL,boundary.value=NULL)
  }
}

check.auglag.blm.active.constraints <- function(object){

  constraints <- object@fit$lambda

  if(any(constraints!=0)){

    i <- which(abs(constraints)!=0)
    U <- object@ineq
    
    list(active.constraints=U[i,],lambda=constraints[i])
  }
  else{
    list(active.constraints=NULL,lambda=NULL)
  }
}

check.auglag.lexpit.active.constraints <- function(object){

  constraints <- object@fit$lambda

  if(any(constraints!=0)){

    i <- which(abs(constraints)!=0)
    U <- cbind(object@ineq$A,object@ineq$B)
    
    list(active.constraints=U[i,],lambda=constraints[i])
  }
  else{
    list(active.constraints=NULL,lambda=NULL)
  }
}

check.lexpit.active.constraints <- function(object,tol=1/1000){

  constraints <- object@constraints$ineq(coef(object))

  if(any(abs(constraints)<tol)){

    i <- which(abs(constraints)<tol)
    U <- cbind(object@ineq$A,object@ineq$B)
    
    list(active.constraints=U[i,],boundary.value=constraints[i])
   }
  else{
    list(active.constraints=NULL,boundary.value=NULL)
  }
}
