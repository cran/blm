\name{logistic.rr}
\alias{logistic.rr}


\title{Estimate a relative risk from a logistic regression model}

\description{

Estimate of relative risk and confidence interval against the reference population of the logistic regression model. This is simply the odds ratio, which approximates the relative risk when the outcome of study is rare.

}

\usage{
logistic.rr(C,fit,alpha=.05,sig=4)
}

\arguments{
  \item{C}{vector of covariate values for the group being compared to the reference population}
   \item{fit}{glm object for logistic regression}
 \item{alpha}{numeric value to specify the 100(1-alpha)\% confidence interval for the risk difference}
  \item{sig}{numeric number of significant digits for reported estimate}
}

\value{
For the relative risk of the group against the reference population, returns a list of \code{est}, \code{lower}, \code{upper} and \code{CI}.
}

\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}
