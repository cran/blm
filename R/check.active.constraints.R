check.active.constraints <- function(object,tol=1/1000){

  constraints <- object@constraints$ineq(coef(object))

  if(any(abs(constraints)<tol)){

    i <- which(abs(constraints)<tol)
    A <- object@constraints$ineq.jac(coef(object))
    
    list(active.constraints=object@ineq[i,],boundary.value=constraints[i],
    			active.grad=A[i,])
  }
  else{
    list(active.constraints=NULL,boundary.value=NULL,active.grad=NULL)
  }
}
