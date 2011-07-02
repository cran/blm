\name{blm.rr}
\alias{blm.rr}


\title{Estimate a relative risk from a binomial linear regression model}

\description{

Estimate of relative risk and confidence interval between subjects of \code{type1} and \code{type2} based on the fit of a binomial linear regression model. 

}

\usage{
blm.rr(type1,type2,fit,sig=4,alpha=.05)
}

\arguments{
  \item{type1}{vector of covariate values for first probability estimate}
  \item{type2}{vector of covariate value for second probability estmate}
   \item{fit}{blm object for linear regression of binary data}
  \item{sig}{numeric number of significant digits for reported estimate}
 \item{alpha}{numeric value to specify the 100(1-alpha)\% confidence interval for the risk difference}
}

\value{
For the relative risk of \code{type1} subject to \code{type2} subject, returns a list of \code{est}, \code{se}, \code{lower}, \code{upper} and \code{CI}.
}

\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\seealso{\link{blm}}