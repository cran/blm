excess.risk <- function(object, group, linear.only=FALSE){
	
	# FUNCTION TO COMPUTE BINNED EXCESS RISK BY GROUPING VARIABLE
	
	if(class(object)[1]!="lexpit"&class(object)[1]!="blm")
		stop("Object must be of class blm or lexpit.")
	
	r <- resid(object)
	
	if(linear.only){
	
		if(class(object)[1]!="lexpit")
			stop("Object must be of class lexpit.")
		
   	 	Z <- model.matrix(object@formula.expit,object@data)
   	 	if(!is.matrix(Z)) Z <- cbind(Z)    
   	   
   	    r <- r+expit(Z%*%object@coef.expit)	
	}
	
	w.bar <- tapply(object@weights, group, sum)
	r.bar <- tapply(object@weights*r, group, sum)
	
r.bar/w.bar
}