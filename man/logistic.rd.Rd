\name{logistic.rd}
\alias{logistic.rd}


\title{Estimate a risk difference between two subject types from a logistic regression model}

\description{

Point estimate and confidence interval for an absolute risk difference based on a logistic regression model of cohort binary data.

}

\usage{
logistic.rd(type1, type2, type1.zero, type2.zero,fit, data, sig = 4, alpha = 0.05)
}

\arguments{
  \item{type1}{vector of covariate values for first probability estimate}
  \item{type2}{vector of covariate value for second probability estmate}
  \item{type1.zero}{vector of logical indicators for type1 covariates to be set to zero}
  \item{type2.zero}{vector of logical indicators for type2 covariates to be set to zero}
  \item{fit}{glm object for logistic regression}
  \item{data}{data frame on which the \code{fit} was based}
  \item{sig}{numeric number of significant digits for reported estimate}
  \item{alpha}{numeric value to specify the 100(1-alpha)\% confidence interval for the risk difference}
}

\value{
For probability of \code{type1} - probability of \code{type2}, returns a list of \code{est}, \code{se}, \code{lower}, \code{upper} and \code{CI}.
}

\details{
For  covariates that are not given a fixed value of set to zero, the empirical distribution is determined from \code{data} and the marginal risk difference is computed by integrating over these factors.
}

\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}
