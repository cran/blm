logistic.rd <- function (type1, type2, fit, data, sig = 4, fixed = FALSE, alpha = 0.05) 
{
    z = qnorm(1 - alpha/2)
    dot.expit = function(x) {
        exp(x)/(1 + exp(x))^2
    }
    if (!fixed) {
        X1 <- model.matrix(fit$formula, data)
        X2 <- X1
        i <- which(type1 != 0)
        j <- which(type2 != 0)
        for (k in i) {
            X1[, k] <- type1[k]
        }
        for (k in j) {
            X2[, k] <- type2[k]
        }
    }
    else {
        X1 <- matrix(type1, nrow = 1)
        X2 <- matrix(type2, nrow = 1)
    }
    pi.1 <- expit(X1 %*% coef(fit))
    pi.2 <- expit(X2 %*% coef(fit))
    var1 <- X1[1, ] %*% vcov(fit) %*% X1[1, ] * dot.expit(X1[1, 
        ] %*% coef(fit))^2
    var2 <- X2[1, ] %*% vcov(fit) %*% X2[1, ] * dot.expit(X2[1, 
        ] %*% coef(fit))^2
    se = sqrt(sum(var1 + var2)/nrow(X1))
    rd = mean(pi.1 - pi.2)
    lower = rd - z * se
    upper = rd + z * se
    CI = paste(round(rd, sig), ", (", round(lower, sig), ", ", 
        round(upper, sig), ")", sep = "", coll = "")
    list(est = rd, se = se, lower = lower, upper = upper, CI = CI)
}


logistic.rr <- function(C,fit,alpha=.05,sig=4){

		est = C%*%coef(fit)
		se = sqrt(t(C)%*%vcov(fit)%*%C)
		z = qnorm(1-alpha/2)
 
    lower = exp(est-z*se)
    upper = exp(est+z*se)
		
		CI = paste(round(exp(est),sig),", (",round(lower,sig),", ",round(upper,sig),")",sep="",coll="")

		list(est=est,se=se,lower=lower,upper=upper,CI=CI)
}


blm.rr <- function(type1,type2,fit,sig=4,alpha=.05){

   z = qnorm(1-alpha/2)
   
    X1 <- matrix(type1,nrow=1)
		X2 <- matrix(type2,nrow=1)  

  pi.1 <- X1%*%coef(fit)
  pi.2 <- X2%*%coef(fit)

  var1 <- X1[1,]%*%vcov(fit)%*%X1[1,]
  var2 <- X2[1,]%*%vcov(fit)%*%X2[1,]

	#VARIANCE FOR LOG TRANFORM
	
  se = sqrt(sum(var1/pi.1^2+var2/pi.2^2))
  est = pi.1/pi.2
  
  lower = exp(log(est)-z*se)
  upper = exp(log(est)+z*se)
  
  CI = paste(round(est,sig),", (",round(lower,sig),", ",round(upper,sig),")",sep="",coll="")

	list(est=est,se=se,lower=lower,upper=upper,CI=CI)
  
}
