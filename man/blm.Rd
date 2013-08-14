\name{blm}
\alias{blm}


\title{Fit a binomial linear regression model}

\description{

A direct probability model for regression with a binary outcome from observational data. 

}

\usage{
blm(formula, data, na.action = na.omit, weights = NULL, 
    strata = NULL, par.init = NULL, warn=FALSE,...)
}

\arguments{
  \item{formula}{formula for linear model for binary outcome, \code{event~x1+x2+...}}
  \item{data}{data.frame containing the variables of \code{formula}}
\item{na.action}{function specifying how missing data should be handled, \link{na.action}}
    \item{weights}{Vector of weights equal to the number of observations. For population-based case-control study, weights are the inverse sampling fractions for controls.}
\item{strata}{vector indicating the stratification for weighted regression with stratified observational data}
\item{par.init}{vector (optional) of initial parameters}
\item{warn}{logical indicator whether to include warnings during algorithm fitting. Default of \code{FALSE} suppresses warnings when testing for feasible parameters.}
\item{...}{Additional arguments passed to \code{constrOptim}}
}

\value{
Returns an object of class \code{blm}.
}

\details{
The \code{blm} model coefficients are the solutions to the maximum of a pseudo log-likelihood using a constrained optimization algorithm with an adaptive barrier method, \code{constrOptim} (Lange, 2010). Variance estimates are based on Taylor linearization (Shah, 2002). When \code{weights} are not NULL, it is assumed that the study is a case-control design.

}

\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\seealso{\link{constrOptim}}

\references{

Kovalchik S, Varadhan R (2013). Fitting Additive Binomial Regression Models with the R Package blm. \emph{Journal of Statistical Software}, 54(1), 1-18. URL: \url{http://www.jstatsoft.org/v54/i01/}.

Lange, K. (2010) \emph{Numerical Analysis for Statisticians}, Springer.

Shah, BV. (2002) Calculus of Taylor deviations. Joint Statistical Meetings.
}

\examples{

data(ccdata)

fit <- blm(y~female+packyear, weights = ccdata$w,strata=ccdata$strata,
				data=ccdata)

summary(fit)


data(aarp)

# ABSOLUTE RISK OF BLADDER CANCER BY 70 YEARS
# FOR DIFFERENT GENDER AND RISK GROUP

fit <- blm(bladder70~female * smoke_status, 
			      data = aarp, 
			      weight=aarp$w)

logLik(fit)

# INTERCEPT IS BASELINE RISK
# ALL OTHER COEFFICIENTS ARE RISK DIFFERENCES FROM BASELINE

summary(fit)

# RISK DIFFERENCE CONFIDENCE INTERVALS (PER 1,000 PERSONS)
confint(fit)*1000
}
