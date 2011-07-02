\name{logistic.rd}
\alias{logistic.rd}


\title{Estimate a risk difference between two subject types from a logistic regression model}

\description{

Point estimate and confidence interval for an absolute risk difference based on a logistic regression model of cohort binary data.

}

\usage{
logistic.rd(type1,type2,fit,data,sig=4,fixed=FALSE,alpha=.05)
}

\arguments{
  \item{type1}{vector of covariate values for first probability estimate}
  \item{type2}{vector of covariate value for second probability estmate}
  \item{fit}{glm object for logistic regression}
  \item{data}{data frame on which the \code{fit} was based}
  \item{sig}{numeric number of significant digits for reported estimate}
  \item{fixed}{logical whether zero-valued quantities should be fixed}
  \item{alpha}{numeric value to specify the 100(1-alpha)\% confidence interval for the risk difference}
}

\value{
For probability of \code{type1} - probability of \code{type2}, returns a list of \code{est}, \code{se}, \code{lower}, \code{upper} and \code{CI}.
}

\details{
If \code{fixed} is \code{FALSE} then the average over all covariates represented as zero in \code{type1} and \code{type2} is taken. If \code{TRUE}, these are fixed at zero and the risk difference applies only to the risk difference between individual who share a zero covriate pattern for these factors. Because of the nonlinearity of the logistic model on the probability scale, the risk difference is a function of all covariates of the model, so will change for different covariate patterns that are still equivalent between the types being compared.
}

\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}
