\name{lexpit}
\alias{lexpit}


\title{Fit a binomial linear-expit regression model}

\description{

A direct probability model for regression with a binary outcome from observational data. Covariate effects are the sum of additive terms and an expit term, which allows for some vectors to be additive and others non-linear.

}

\usage{
lexpit(f.linear,f.expit,data,par.init,weights=NULL,ineq=NULL,trace=FALSE,
	tol=1e-6,augmented=TRUE,warn=-1,...)
}

\arguments{
  \item{f.linear}{formula for linear model for binary outcome, \code{event~x1+x2+...}}
  \item{f.expit}{formula for expit model, linear in expit, \code{event~z1+z2+...}}
  \item{data}{data.frame containing the variables of \code{f.linear} and \code{f.expit}}
  \item{par.init}{List of \code{linear} and \code{expit} logical or numeric.}
    \item{weights}{Vector of weights equal to the number of observations. Weights can be used to estimate population statistics from survey or case-control data.}
  \item{ineq}{optional user-specified inequality constraint matrix.}
\item{trace}{logical, whether optimization details at each outer iteration should be printed.}
\item{tol}{numeric criterion for judging a constraint active when \code{augmented} is \code{FALSE}.}
\item{augmented}{logical, whether the optimization method should use constrained optimization with an augmented Lagrangian (see details).}
\item{warn}{numerical, passed to \code{options} during run. Default suppresses warning messages.}
\item{...}{Additional argument to \code{control.outer} of \code{constrOptim.nl}}
}

\value{
Returns an object of class \code{lexpit}.
}

\details{
\code{lexpit} like \code{blm} is a wrapper for \code{auglag} and \code{constrOptim.nl} from the \code{alabama} package, which performs constrained optimization for nonlinear functions and linear/nonlinear inequality constraints. \code{constrOptim.nl} performs the constrained optimization procedure with a interior point adaptive barrier (Lange, 2010). If augmented is TRUE than an augmented Lagrangian method is used (Madsen, 2004).


When \code{par.init} is not specified an iterative search for feasible starting values using a linear projection from the OLS estimates into the feasible region. If a list is given, it has the named components \code{linear} and \code{expit} where each is a vector of the same length as the linear and expit parameters. The initial parameters are set equal to the supplied vectors if these are numeric. If logical, they are taken as instructions for setting the specified parameter to zero when \code{TRUE} and to use the iterative procedure when \code{FALSE}. 

The \code{ineq} matrix for the linear component of the model is determined automatically from the unique covariate patterns of the design matrix contains only the extreme values (endpoints of empirical range) for each column. It can also be supplied directly by the user. In this case, it is a matrix with the rows equal to the number of inequality constraints for the lower bound and the columns corresponding to each column of the model design matrix.

When a constraint is at the boundary of the feasible region, the constraint will be approximately zero when evaluated at the final estimates of the algorithm. When using the augmented Lagrangian, an active constraint is defined as one where the multiplier \code{lambda} is non-zero. When \code{augmented} is \code{FALSE} (i.e. classical barrier method is used) \code{tol} is the criterion for zero, so any evaluated constraints less than or equal to \code{tol} are considered active. If any constraints are active, the covariance-variance based on the Taylor series approximation might be inaccurate. Since the augmented Lagrangian procedure takes into the active constraints into account when computing the Hessian matrix, whose inverse is the covariance-variance of the model estimates, this is the recommended method for estimation at the boundary.

}

\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\seealso{\link{blm-class},\link{blm},\link{summary},\link{constrOptim.nl}}

\references{
Lange, K (2010). \emph{Numerical Analysis for Statisticians}, Springer.

Madsen K, Nielsen HB, Tingleff O (2004). \emph{Optimization with Constraints}, IMM, DTU,
}

\examples{

data(grad)

# EFFECT OF GRE ON UNDERGRAD SCHOOL PRESTIGE ON GRAD ADMISSION
# CONTROLLED FOR SCHOOL RANK

fit <- lexpit(admit~I(scale(gre)),admit~factor(rank),grad)

fit 

summary(fit)

# PROB ADMISSION FOR GRE 2 SD ABOVE AVERAGE
# ADJUSTED FOR MOST PRESTIGIOUS UNDERGRAD SCHOOL

is.constrained(fit)

ci(fit,C=list(linear=2,expit=c(1,0,0,0)),coef=FALSE,
	average=list(linear=FALSE,expit=c(FALSE,FALSE,FALSE,FALSE)),
	method="logit") 


### USER-SUPPLIED INITIAL VALUES

fit2 <- lexpit(admit~I(scale(gpa)),admit~factor(rank),grad,
			par.init=list(linear=0,expit=runif(4)))

fit2

### IF FIT WITH ADAPTIVE BARRIER METHOD
### STANDARD ERROR MIGHT BE INCORRECT WHEN CONSTRAINTS ARE ACTIVE

fit2 <- lexpit(admit~I(scale(gpa)),admit~factor(rank),grad,augmented=FALSE)

fit2

}
