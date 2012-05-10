# POP CASE-CONTROL SIM
library(blm)

n <- 500

data <- data.frame(
	female = rep(c(0,1),length=n),
	packyear = sample(c(0,10,20,30),n,rep=TRUE,prob=c(.3,.1,.3,.2)),
	strata = rep(1:5,length=n)
)

data$y <- sapply(data$female*.05+expit(logit(.15)+.05*data$packyear),
					function(x)rbinom(1,1,x))

ncases <- tapply(data$y,data$strata,sum)
nstrata <- tapply(data$y,data$strata,length)

# CREATE CASE-CONTROL DATASET
index <- list()

for(i in 1:length(ncases)){
		index[[i]] <- sample(row.names(data)[data$strata==i&data$y==0],ncases[i])		
}

data <- data[c(unlist(index),which(data$y==1)),]
data$w <- (nstrata/ncases)[factor(data$strata)]
data$w[data$y==1] <- 1

formula.linear <- y~female
formula.expit <- y~1
weights <- data$w
strata <- data$strata

fit <- lexpit(formula.linear,formula.expit,data=data,weights=weights,strata=strata)
