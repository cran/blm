gof <- function(object){
	
	if(class(object)[1]=="lexpit"|class(object)[1]=="blm"){
		p <- predict(object)
		y <- object@y
		w <- object@weights
	}
	else{
		p <- expit(predict(object)) # LOGISTIC
		y <- object$y
		w <- object$prior.weights
	}	
		brks <- quantile(p, seq(0, 1, length=11))
		g <- cut(p, br = brks, include.lowest = TRUE)
		O <- tapply(y,g,sum)
		E <- tapply(p*w,g,sum)

		result <- data.frame(O=O,E=E)
	    chisq = sum((O-E)^2/E)
	    P = 1 - pchisq(chisq, df=8)


list(table = result, chisq = chisq, p.value = P)
}

EO <- function(object, index=NULL,level=.95){
	if(class(object)[1]!="lexpit"&class(object)[1]!="blm")
		stop("Object must be an instance of a blm or lexpit model.")
	
	if(is.null(index)) index <- rep("Overall", length(object@y))
	
	p <- predict(object)
	
	O <- tapply(object@y, index, sum)
	E <- tapply(p*object@weights, index, sum)
	
	 z <- qnorm(1-(1-level)/2)
	 lower <- E/O*exp(-z*sqrt(1/O))
     upper <- E/O*exp(+z*sqrt(1/O))
      
  	results <- data.frame(
  					E = E, 
  					O = O, 
  					EtoO = E/O,
  				    lowerCI = lower, 
  				    upperCI = upper)
 
 results
}