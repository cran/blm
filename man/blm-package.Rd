\name{blm-package}
\alias{blm-package}
\docType{package}
\title{
Binomial linear and linear-expit regression  model
}
\description{
The functions \code{blm} and \code{lexpit} implement a binomial linear and linear-expit regression model. Estimates are the maximum likelihood estimates with constrained optimization through an interior-point method or an Augmented Lagrangian Adaptive Barrier Minimization Algorithm to ensure that estimable probabilities are in the (0,1) interval. 
}
\details{
\tabular{ll}{
Package: \tab blm\cr
Type: \tab Package\cr
Version: \tab 2011.1.1\cr
Date: \tab 2011-08-15\cr
Depends: R (>= 2.10.1), stats, BB, MASS, numDeriv, alabama (>= 2010.9-1)\cr
License: \tab GPL (>= 2) \cr
LazyLoad: \tab yes
}

The \code{blm} package extends the \code{alabama} package of Ravi Varadhan.
}

\author{

Maintainer: S. Kovalchik <s.a.kovalchik@gmail.com>

}

\seealso{\link{alabama}}

\keyword{ package }
