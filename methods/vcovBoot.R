boot.blm <- function(object,n.boot=25){

		mle = object@fit$par
		p = length(mle)
		boot.estimates <- matrix(0,n.boot,p)
		N <- nrow(object@data)

                i = 0      
		while(i < n.boot){

                  sample <- sample(1:N,replace=TRUE)
		  data <- object@data[sample,]
                                
		  one.step.fit <-  tryCatch(blm(object@formula,data,	   	                                                   ineq = object@ineq,
                                par.init=mle,
                                control.outer=list(trace=FALSE)),
                                error=function(e){NA})
                  
                  if(class(one.step.fit)!="logical"){
                    i = i + 1
                    boot.estimates[i,] <- one.step.fit@fit$par
                   }
		}
                                
   boot.estimates
}

boot.lexpit <- function(object,n.boot=25){

				mle = object@fit$par
                                fit = summary(object)
                                gamma = fit$est.expit
                                beta = fit$est.linear
				p = length(mle)
				boot.estimates <- matrix(0,n.boot,p)
				N <- nrow(object@data)

                i = 0      
		while(i < n.boot){
		  sample <- sample(1:N,replace=TRUE)
		  data <- object@data[sample,]
                                
		  one.step.fit <-  tryCatch(lexpit(object@formula.linear,object@formula.expit,data,	   	                          ineq = object@ineq,
                                       par.init=list(linear=beta,expit=gamma),
                                       control.outer=list(trace=FALSE)),
                                            error=function(e){NA})
                  
                  if(class(one.step.fit)!="logical"){
                    i = i + 1
                    boot.estimates[i,] <- one.step.fit@fit$par
                   }
		}
				
  boot.estimates
}

subboot.blm <- function(object,n.boot=25,alpha=runif(1,.5,1)){

		mle = object@fit$par
		p = length(mle)
		boot.estimates <- matrix(0,n.boot,p)
		N <- nrow(object@data)

                i = 0      
		while(i < n.boot){

                  sample <- sample(1:N,size=floor(N^alpha))
		  data <- object@data[sample,]
                                
		  one.step.fit <-  tryCatch(blm(object@formula,data,	   	                                                   ineq = object@ineq,
                                par.init=mle,
                                control.outer=list(trace=FALSE)),
                                error=function(e){NA})
                  
                  if(class(one.step.fit)!="logical"){
                    i = i + 1
                    boot.estimates[i,] <- one.step.fit@fit$par
                   }
		}
                                
   boot.estimates
}

subboot.lexpit <- function(object,n.boot=25,alpha=runif(1,.5,1)){

				mle = object@fit$par
                                fit = summary(object)
                                gamma = fit$est.expit
                                beta = fit$est.linear
				p = length(mle)
				boot.estimates <- matrix(0,n.boot,p)
				N <- nrow(object@data)

                i = 0      
		while(i < n.boot){
		  sample <-  sample(1:N,size=floor(N^alpha))
		  data <- object@data[sample,]
                                
		  one.step.fit <-  tryCatch(lexpit(object@formula.linear,object@formula.expit,data,	   	                          ineq = object@ineq,
                                       par.init=list(linear=beta,expit=gamma),
                                       control.outer=list(trace=FALSE)),
                                            error=function(e){NA})
                  
                  if(class(one.step.fit)!="logical"){
                    i = i + 1
                    boot.estimates[i,] <- one.step.fit@fit$par
                   }
		}
				
  boot.estimates
}



setGeneric("boot",function(object,n.boot=25){standardGeneric("boot")})

setMethod("boot","blm",boot.blm)
setMethod("boot","lexpit",boot.lexpit)

setGeneric("subboot",function(object,n.boot=25,alpha=runif(1,.5,1)){standardGeneric("subboot")})

setMethod("subboot","blm",subboot.blm)
setMethod("subboot","lexpit",subboot.lexpit)
