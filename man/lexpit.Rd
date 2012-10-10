\name{lexpit}
\alias{lexpit}


\title{Fit a linear-expit regression model}

\description{

A direct probability model for regression with a binary outcome from observational data. Covariate effects are the sum of additive terms and an expit term, which allows some explanatory variables to be additive and others non-linear.

}

\usage{
lexpit(formula.linear,formula.expit,data,na.action=na.omit,
					weights=NULL,strata=NULL,par.init=NULL,
					warn = FALSE,
					control.lexpit=list(max.iter=1000,tol=1E-7),...)
}

\arguments{
  \item{formula.linear}{formula for linear model for binary outcome, \code{event~x1+x2+...}}
  \item{formula.expit}{formula for expit model, linear in expit, \code{event~z1+z2+...}}
  \item{data}{data.frame containing the variables of \code{formula.linear} and \code{formula.expit}}
\item{na.action}{function specifying how missing data should be handled, \link{na.action}}
    \item{weights}{Vector of weights equal to the number of observations. For population-based case-control study, weights are the inverse sampling fractions for controls.}
\item{strata}{vector indicating the stratification for weighted regression with stratified observational data}
\item{par.init}{list (optional) of initial parameters for \code{linear} and \code{expit} terms.}
\item{warn}{logical indicator whether to include warnings during algorithm fitting. Default of \code{FALSE} suppresses warnings when testing for feasible parameters.}
  \item{control.lexpit}{list with control parameters for optimization algorithm}
\item{...}{Additional arguments passed to \code{constrOptim}}
}

\value{
Returns an object of class \code{lexpit}.
}

\details{
\code{lexpit} model uses a two-stage optimization procedure. At the first stage linear terms the solutions to the maximum of a pseudo log-likelihood using a constrained optimization algorithm with an adaptive barrier method, \code{constrOptim} (Lange, 2010). The second stage maximizes the pseudo log-likelihood with respect to the expit terms using iterative reweighted least squares with an offset term for the linear component of the model.

Variance estimates are based on Taylor linearization (Shah, 2002). When \code{weights} are not NULL, it is assumed that the study is a case-control design. 

}

\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\seealso{\link{constrOptim}, \link{nlm}}

\references{

Lange, K. (2010) \emph{Numerical Analysis for Statisticians}, Springer.


Shah, BV. (2002) Calculus of Taylor deviations. Joint Statistical Meetings.
}

\examples{

data(ccdata)

fit <- lexpit(y~female,y~packyear,weights = ccdata$w,
       		strata=ccdata$strata,data=ccdata)

summary(fit)

# LEXPIT MODEL FOR BLADDER CANCER RISK BY AGE 70
formula.linear <- bladder70~female * smoke_status
formula.expit <- bladder70~redmeat+fiber.centered+I(fiber.centered^2)

# ADDITIVE EFFECTS FOR GENDER AND SMOKING
# LOGISTIC EFFECTS FOR FIBER AND REDMEAT CONSUMPTION
data(aarp)

fit <- lexpit(formula.linear, formula.expit, aarp, weight=aarp$w)

# SUMMARY
summary(fit)
confint(fit)

# FITTED ABSOLUTE RISK PER 1,000 PERSONS
head(predict(fit)*1000)
}
