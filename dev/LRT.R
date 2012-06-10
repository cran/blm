# HYPOTHESIS TESTING FUNCTIONS FOR BLM
LRT <- function(object,...){
	
	f <- object@formula
	
	labels <- attr(terms(f),"term.labels")
	formulas <- lapply(labels, function(term) update(f, paste("~.-",term, sep="")))
	
	LL <- sapply(formulas, function(new.formula) blm(new.formula, object@data)@loglik)
	
	LRTs <- 2*(object@loglik-LL)
	
	table <- cbind(LRT=LRTs, pvalue=1-pchisq(LRTs,df=1))
	row.names(table) <- labels

table
}

parametric_bootstrap <- function(fit.big, fit.small, n.boot=1000){
	
	# SIMULATE LRT WITH PARAMETRIC BOOTSTRAP
	p <- predict(fit.small) # PARAMETRIC FORM
	
	f.big <- update(fit.big@formula, simy~.)
	f.small <- update(fit.small@formula, simy~.)

	obs.T <- 2*(fit.big@loglik-fit.small@loglik)	
	
	sim.T <- sapply(rep(1,n.boot), function(void){
		
		simy <- rbinom(length(p), size=1, prob=p)
		
		fit.big@data$simy <- simy
		fit.small@data$simy <- simy
		
		fit.sim.big <- blm(f.big, fit.big@data)
		fit.sim.small <- blm(f.small, fit.small@data)
		
	 2*(fit.sim.big@loglik-fit.sim.small@loglik)	
	})

sim.T	
}

# TESTING
library(blm)
library(survival)

data(mgus)

# BLM
mgus$event <- ifelse(mgus$futime<5000&mgus$death==1,1,0)
mgus <- mgus[!is.na(mgus$mspike),]

fit.big <- blm(event~I(scale(age))+mspike+sex, data=mgus)
fit.small <- blm(event~I(scale(age))+sex, data=mgus)

obs.T <- 2*(fit.big@loglik-fit.small@loglik)

sim.T <- parametric_bootstrap(fit.big, fit.small, 1000)

(sum(sim.T>=obs.T)+1)/(length(sim.T)+1)


# LEXPIT
fit.big <- lexpit(event~sex,event~I(scale(age)), data=mgus)
fit.small <- blm(event~I(scale(age))+sex, data=mgus)

obs.T <- 2*(fit.big@loglik-fit.small@loglik)

sim.T <- parametric_bootstrap(fit.big, fit.small, 1000)

(sum(sim.T>=obs.T)+1)/(length(sim.T)+1)



# SIMULATE A DATASET TO CHECK LARGE-SAMPLE AGREEMENT
test <- data.frame(
				x = runif(1000),
				noise = runif(1000)
)

test$event <- rbinom(n=nrow(test),size=1,prob=0.1+test$x*.5)

fit <- blm(event~x+noise, data=test)


sim <-   sapply(rep(1,500), function(void) {
					
		test <- data.frame(
				x = runif(1000),
				noise = runif(1000)
			)

		test$event <- rbinom(n=nrow(test),size=1,prob=0.1+test$x*.5)

		fit <- blm(event~x+noise, data=test)

		fit@coef
})

beta.bar <- apply(sim,1,mean)
diff <- sim-matrix(beta.bar,nr=dim(sim)[1],nc=dim(sim)[2])
V <- apply(diff, 2, function(x) outer(x, x))
V <- matrix(rowSums(V),length(beta.bar), length(beta.bar))/ncol(sim)

# JACKKNIFE SIMULATION
betas <-  sapply(1:nrow(test), function(i) blm(event~x+noise, data=test[-i,])@coef)
n <- ncol(betas)
pseudo <- matrix(n*coef(fit),nrow=3,ncol=ncol(betas))-(n-1)*betas
diff.betas <- pseudo-matrix(apply(pseudo,1,mean), nr=3, nc=ncol(betas))
V <- apply(diff.betas, 2, function(x) outer(x,x))
V <- matrix(rowSums(V), nrow(betas), nrow(betas))/(ncol(betas)*(ncol(betas)-1))

V <- rowSums(diff.betas*diff.betas)/(ncol(betas)*(ncol(betas)-1))
sqrt(diag(vcov(fit)))

# SIMULATION TO CHECK LARGE-SAMPLE AGREEMENT
test <- data.frame(
				x1 = runif(500),
				x2 = runif(500),
				noise = runif(500)
)

test$event <- rbinom(n=nrow(test),size=1,prob=0.1*test$x1+expit(logit(0.25)+test$x2))
fit <- lexpit(event~x1, event~x2, data=test)
summary(fit)
LRT(fit)

fit.logit <- glm(event~x1+x2,fam="binomial",data=test)

sim <-   sapply(rep(1,500), function(void) {
					
		test <- data.frame(
				x1 = runif(500),
				x2 = runif(500),
				noise = runif(500)
		)

		test$event <- rbinom(n=nrow(test),size=1,prob=0.1*test$x1+expit(logit(0.25)+test$x2))
		fit <- lexpit(event~x1, event~x2, data=test)

		c(fit@coef.linear,  fit@coef.expit)
				})

beta.bar <- apply(sim,1,mean)
diff <- sim-matrix(beta.bar,nr=dim(sim)[1],nc=dim(sim)[2])
V <- apply(diff, 2, function(x) outer(x, x))
V <- matrix(rowSums(V),length(beta.bar), length(beta.bar))/ncol(sim)
