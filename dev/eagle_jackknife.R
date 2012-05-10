### BOOTSTRAP FIT
eagle.sample <- function(index,f.linear,f.expit,data){
	
	stratum <- data$strata[index]
	data <- data[-index,]
	data$sf_1[data$strata==stratum] <- (data$nstrata/(data$nstrata-1))[data$strata==stratum]*data$sf_1[data$strata==stratum]
	
	fit <- lexpit(f.linear,f.expit,weights=data$sf_1,strata=NULL,data=data)

coef(fit)
}