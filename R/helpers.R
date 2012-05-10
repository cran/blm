expit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))
dexpit <- function(x) expit(x)*(1-expit(x))
ddexpit <- function(x) dexpit(x)*(1-exp(x))-expit(x)*exp(x)

