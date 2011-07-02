blm.ineq.matrix <- function(f,data){

 is.continuous <- function(X){
  apply(X,2,function(x){length(unique(x))>2})
 }

 any.interaction <- function(X){
  names <- colnames(X)
  length(grep(":",names))>0
 }

 which.interaction <- function(X){
  names <- colnames(X)
  which.interactions <- grep(":",names)
   interactions <- lapply(which.interactions,function(i){
      strsplit(names[i],":")[[1]]
    })

  names(interactions) <- names[which.interactions]

  interactions
 }

##############
 
       X <- model.matrix(f,data)
       if(!is.matrix(X)) X = matrix(X,ncol=1)
       var.continuous <- is.continuous(X)

       if(any(var.continuous)){
         j <- which(var.continuous)
         
       for(i in j){
           X[,i] <- sample(rep(range(X[,i]),length=nrow(X)))
          }
       }

       if(any.interaction(X)){
         interaction.list <- which.interaction(X)
           for(i in names(interaction.list)){
            
             X[,i] <- X[,interaction.list[[i]][1]]*X[,interaction.list[[i]][2]]
           }
       }
       
   unique(X)
}


blm.constraints <- function(f,data,A=NULL){
  
  if(is.null(A)){
    A <- blm.ineq.matrix(f,data)
  }
  
  ineq <- function(par){c(A%*%par,1-A%*%par)}

  ineq.jac <- function(par){
    
     rbind(A,-A)
  }

  return(list(A=A,ineq=ineq,ineq.jac=ineq.jac))
}


lexpit.ineq.matrix <- function(f,data){
		
        A <- blm.ineq.matrix(f,data)
        
        if(attr(terms(f),"int")) A = A[,-1]
        if(!is.matrix(A)) A = matrix(A,ncol=1)
   A
}


lexpit.constraints <- function(f.linear,f.expit,data,A=NULL){
  
  if(is.null(A)){
    A <- lexpit.ineq.matrix(f.linear,data)
  }
  
  X <- model.matrix(f.linear,data)
  Z <- model.matrix(f.expit,data)

  if(attr(terms(f.linear),"int")) X = X[,-1]
  p = ifelse(is.matrix(X),ncol(X),1)
  X = matrix(X,ncol=p)
  q = ncol(Z)
  
  ineq <- function(par){
    
    beta = par[1:p]
    gamma = par[(p+1):(p+q)]
    Z.range <- range(expit(Z%*%gamma))
    
    c(A%*%beta+Z.range[1],1-A%*%beta-Z.range[2])
  }

  ineq.jac <- function(par){
    
     beta = par[1:p]
     gamma = par[(p+1):(p+q)]
     
     pi.gamma <- expit(Z%*%gamma)
     min = min(pi.gamma)
     max = max(pi.gamma)
     
     i.min = which(pi.gamma==min)[1]
     i.max = which(pi.gamma==max)[1]
     
     Z.min = matrix(0,nrow(A),q)
     Z.max = Z.min
     
     dot.expit = function(x){exp(x)/(1+exp(x))^2}

     for(i in 1:nrow(A)){
     	Z.min[i,] <- Z[i.min,]*dot.expit(min)
     	Z.max[i,] <- Z[i.max,]*dot.expit(max)
     }
     
     gamma.jac <- rbind(Z.min,-Z.max)
     beta.jac = rbind(A,-A)
     cbind(beta.jac,gamma.jac)    
  }

  return(list(A=A,ineq=ineq,ineq.jac=ineq.jac))

}

