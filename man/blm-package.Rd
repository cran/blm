\name{blm-package}
\alias{blm-package}
\docType{package}
\title{
Binomial linear and linear-expit regression  model
}
\description{
The functions \code{blm} and \code{lexpit} implement a binomial linear and linear-expit regression model. Estimates are the maximum likelihood estimates with constrained optimization through adaptive barrier method to ensure that estimable probabilities are in the (0,1) interval. 
}
\details{
\tabular{ll}{
Package: \tab blm\cr
Type: \tab Package\cr
Version: \tab 2013.2.4.4\cr
Date: \tab 2013-8-14\cr
Depends: \tab R (>= 2.10.1), methods\cr
Imports: \tab stats, stats4\cr
License: \tab GPL (>= 2) \cr
LazyLoad: \tab yes
}

}

\author{

Maintainer: Stephanie Kovalchik <s.a.kovalchik@gmail.com>

}

\references{

Kovalchik S, Varadhan R (2013). Fitting Additive Binomial Regression Models with the R Package blm. \emph{Journal of Statistical Software}, 54(1), 1-18. URL: \url{https://www.jstatsoft.org/v54/i01/}.

}

\seealso{\link{constrOptim}, \link{blm}, \link{lexpit}}

\keyword{ package }
