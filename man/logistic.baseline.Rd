\name{logistic.baseline}
\alias{logistic.baseline}


\title{Estimate an absolute risk from a logistic regression model}

\description{

Estimate the absolute risk and confidence interval against the reference population of the logistic regression model of cohort data.

}

\usage{
logistic.baseline(C,fit,sig=4,alpha=.05)
}

\arguments{
 \item{C}{vector of covariate values for the group whose absolute risk is to be estimated}
 \item{fit}{glm object for logistic regression}
 \item{sig}{numeric number of significant digits for reported estimate}
 \item{alpha}{numeric value to specify the 100(1-alpha)\% confidence interval for the risk difference}

}

\value{
For the relative risk of the group against the reference population, returns a list of \code{est}, \code{lower}, \code{upper} and \code{CI}.
}

\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}
