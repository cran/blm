\name{blm}
\alias{blm}


\title{Fit a binomial linear regression model}

\description{

Binomial linear model for cohort binary data. 

}

\usage{
blm(f,data,par.init,ineq=NULL,trace=FALSE,tol=1e-6,augmented=TRUE,warn=-1,...)
}

\arguments{
  \item{f}{formula for linear model for binary outcome, \code{event~x1+x2+...}}
  \item{data}{data.frame containing the variables of \code{f}}
  \item{par.init}{optional vector of user-specified initial parameters}
  \item{ineq}{optional user-specified inequality constraint matrix.}
\item{trace}{logical, whether details of optimization should be printed with each outer iteration.}
\item{tol}{numeric criterion for judging a constraint active when \code{augmented} is \code{FALSE}.}
\item{augmented}{logical, whether the optimization method should use constrained optimization with an augmented Lagrangian (see details).}
\item{warn}{numerical, passed to \code{options} during run. Default suppresses warning messages.}
\item{...}{Additional argument to \code{control.outer} of \code{auglag} and \code{constrOptim.nl}}
}

\value{
Returns an object of class \code{blm}.
}

\details{
\code{blm} is a wrapper for \code{auglag} and \code{constrOptim.nl} from the \code{alabama} package, which performs constrained optimization for nonlinear functions and linera/nonlinear inequality constraints. \code{constrOptim.nl} performs the constrained optimization procedure with a interior point adaptive barrier (Lange, 2010). If augmented is TRUE than an augmented Lagrangian method is used (Madsen, 2004).

When \code{par.init} is not specified, a linear projection of the OLS estimators into the feasible region are used.

The \code{ineq} matrix is determined automatically from the unique covariate patterns of the design matrix containing only the extreme values (endpoints of empirical range) for each column of continuous variables. It can also be supplied directly by the user. In this case, it is a matrix with the rows equal to the number of inequality constraints for the lower bound and the columns corresponding to each column of the model design matrix.


When a constraint is at the boundary of the feasible region, the constraint will be approximately zero when evaluated at the final estimates of the algorithm. When using the augmented Lagrangian, an active constraint is defined as one where the multiplier \code{lambda} is non-zero. When \code{augmented} is \code{FALSE} (i.e. classical barrier method is used) \code{tol} is the criterion for zero, so any evaluated constraints less than or equal to \code{tol} are considered active. If any constraints are active, the covariance-variance based on the Taylor series approximation might be inaccurate. Since the augmented Lagrangian procedure takes into the active constraints into account when computing the Hessian matrix, whose inverse is the covariance-variance of the model estimates, this is the recommended method for estimation at the boundary.
}


\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\seealso{\link{blm-class}, \link{constrOptim.nl},\link{summary}}

\references{
Lange, K (2010). \emph{Numerical Analysis for Statisticians}, Springer.

Madsen K, Nielsen HB, Tingleff O (2004). \emph{Optimization with Constraints}, IMM, DTU,
}
\examples{

data(grad)

fit <- blm(admit~factor(rank),grad)

fit 

summary(fit)

ci(fit,c(1,0,0,0))  #PROB GRAD SCHOOL ADMISSION 
		    #FOR STUDENTS FROM MOST PRESTIGIOUS SCHOOL


### INCLUDE FACTORS FOR UNDERGRAD GPA AND GRE AND SUPPLY INITIAL VALUES

fit2 <- blm(admit~I(scale(gre))+I(scale(gpa))+factor(rank),grad)

fit2

summary(fit2)

### IMPROVEMENT IN AIC USING 2 DEGREES OF FREEDOM

summary(fit)$AIC-summary(fit2)$AIC 

### IF FIT WITH ADAPTIVE BARRIER METHOD
### STANDARD ERROR MIGHT BE INCORRECT WHEN CONSTRAINTS ARE ACTIVE

fit2 <- blm(admit~I(scale(gre))+I(scale(gpa))+factor(rank),grad,augmented=FALSE)

fit2

}