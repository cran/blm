dispersion.blm <- function(object){

     if(class(object)=="lexpit"){
  	X <- model.matrix(object@formula.linear,object@data)
        Z <- model.matrix(object@formula.expit,object@data)
         if(attr(terms(object@formula.linear),"intercept")==1){
           X = X[,-1]
         }
         if(!is.matrix(X)){
           X = matrix(X,ncol=1)
         }
        X = cbind(X,Z)
	Y <- model.frame(object@formula,object@data)[,1]
	covariate.class <- apply(X,1,function(x){paste(x,sep="",collapse="")})

        }
     else{
         
	X <- model.matrix(object@formula,object@data)
	Y <- model.frame(object@formula,object@data)[,1]

        }
  
	covariate.class <- apply(X,1,function(x){paste(x,sep="",collapse="")})
	
	pi = predict(object)

	O = tapply(Y,covariate.class,sum)

	N = tapply(Y,covariate.class,length)
	E = N*tapply(pi,covariate.class,mean)
	
	n = length(O)
	p = length(as.vector(coef(object)))
	X2 = sum((O-E)^2/E)
	pearson.df = n-p
	pearson.p = 1-pchisq(X2,df=pearson.df)
	
	D = 2*(sum((Y*log(Y/pi))[Y!=0])+sum(((1-Y)*log((1-Y)/(1-pi)))[Y!=1]))
	deviance.df = sum(N)-p
	deviance.p = 1-pchisq(D,df=deviance.df)

	list(
		observed = O,
		expected = E,
		deviance=D,
		pearson=X2,
		pearson.df=pearson.df,
		deviance.df=deviance.df,
		pearson.p=pearson.p,
		deviance.p=deviance.p
		)

}

setGeneric("dispersion",function(object){standardGeneric("dispersion")})

setMethod("dispersion","blm",dispersion.blm)
setMethod("dispersion","lexpit",dispersion.blm)


logistic.dispersion <- function(f,data,fit){

	X <- model.matrix(f,data)
        if(!is.matrix(X)) X = matrix(X,ncol=1)
    
	Y <- model.frame(f,data)[,1]
        
	covariate.class <- apply(X,1,function(x){paste(x,sep="",collapse="")})
	
	pi = expit(X%*%coef(fit))

	O = tapply(Y,covariate.class,sum)

	N = tapply(Y,covariate.class,length)
	E = N*tapply(pi,covariate.class,mean)
	
	n = length(O)
	p = length(as.vector(coef(fit)))
	X2 = sum((O-E)^2/E)
	pearson.df = n-p
	pearson.p = 1-pchisq(X2,df=pearson.df)
	
	D = 2*(sum((Y*log(Y/pi))[Y!=0])+sum(((1-Y)*log((1-Y)/(1-pi)))[Y!=1]))
	deviance.df = sum(N)-p
	deviance.p = 1-pchisq(D,df=deviance.df)

	list(
		observed = O,
		expected = E,
		deviance=D,
		pearson=X2,
		pearson.df=pearson.df,
		deviance.df=deviance.df,
		pearson.p=pearson.p,
		deviance.p=deviance.p
		)

}
