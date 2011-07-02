vcovBoot.blm <- function(object,n.boot=25){

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
                                control.outer=list(trace=FALSE,itmax=1)),
                                            error=function(e){NA})
                  
                  if(class(one.step.fit)!="logical"){
                    i = i + 1
                    boot.estimates[i,] <- one.step.fit@fit$par
                   }
		}
                                
   matrix(apply(apply(boot.estimates,1,function(x){outer(x-mle,x-mle)}),1,mean),p,p)
}

vcovBoot.lexpit <- function(object,n.boot=25){

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
                                
		  one.step.fit <-  tryCatch(lexpit(object@formula.linear,object@formula.expit,data,	   	                           ineq = object@ineq,
                                par.init=list(linear=beta,expit=gamma),
                                control.outer=list(trace=FALSE,itmax=1)),
                                            error=function(e){NA})
                  
                  if(class(one.step.fit)!="logical"){
                    i = i + 1
                    boot.estimates[i,] <- one.step.fit@fit$par
                   }
		}
				
  matrix(apply(apply(boot.estimates,1,function(x){outer(x-mle,x-mle)}),1,mean),p,p)
}

setGeneric("vcovBoot",function(object,n.boot=25){standardGeneric("vcovBoot")})

setMethod("vcovBoot","blm",vcovBoot.blm)
setMethod("vcovBoot","lexpit",vcovBoot.lexpit)
