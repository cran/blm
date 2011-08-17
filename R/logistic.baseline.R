logistic.baseline <- function(C,fit,sig = 4,alpha = 0.05) {
 
   z = qnorm(1 - alpha/2)

    dot.expit = function(x) {
        exp(x)/(1 + exp(x))^2
    }
	
    beta = coef(fit)
    se = sqrt(C%*%vcov(fit)%*%C*dot.expit(C%*%beta)^2)
    est = expit(C%*%beta)

    lower = est - z * se
    upper = est + z * se

    CI = paste(round(est, sig), ", (", round(lower, sig), ", ", 
        round(upper, sig), ")", sep = "", coll = "")

    list(est = est, se = se, lower = lower, upper = upper, CI = CI)
}
