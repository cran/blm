excess.risk <- function(object, group){
	
	# FUNCTION TO COMPUTE BINNED EXCESS RISK BY GROUPING VARIABLE
	
	if(class(object)[1]!="lexpit"&class(object)[1]!="blm"){
		if(object$family$link!="logit") stop("Object must be of class blm, lexpit, or glm with logit link.")
			excess.logit.risk(object, group)
	}
	else{
	p <- predict(object)
	p.bar <- tapply(p, group, mean)
	w.bar <- tapply(object@weights, group, sum)
	y.bar <- tapply(object@weights*object@y, group, sum)
	
	y.bar/w.bar-p.bar
	}	
}



excess.logit.risk <- function(object, group){
	
	w.bar <- tapply(object$prior.weights, group, sum)
	p.bar <- tapply(expit(predict(object)), group, sum)
	y.bar <- tapply(object$y*object$prior.weights, group, sum)

y.bar/w.bar-p.bar
}


risk.exposure.plot <- function(y, group, weights=NULL,...){
	
risk <- function(start, stop, y, group, weights){
	lower.upper <- quantile(group, c(start, stop))
	index <- which(group>=lower.upper[1]&group<=lower.upper[2])
sum(y[index])/sum(weights[index])
}

covariate <- function(start, stop, group, weights){

	lower.upper <- quantile(group, c(start, stop))
	index <- which(group>=lower.upper[1]&group<=lower.upper[2])

sum(group[index]*weights[index])/sum(weights[index])
}

	if(is.null(weights)) weights <- rep(1, length(y))
	
Y <- mapply(risk, start=seq(0,.8,by=.01), stop=seq(.2,1,by=.01), MoreArgs=list(group=group, y=y, weights=weights))
X <- mapply(covariate, start=seq(0,.8,by=.01), stop=seq(.2,1,by=.01), MoreArgs=list(group=group, weights=weights))


data <- list(
	risk = Y,
	x = X
)

	invisible(data)
	
	scatter.smooth(y=data$risk, x=data$x, ylab="Average risk",...)

}