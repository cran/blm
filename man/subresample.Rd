\name{subresample}
\docType{methods}
\alias{subresample}
\alias{subresample,blm-method}
\alias{subresample,lexpit-method}

\title{Get a resampled subsample of \code{blm} and \code{lexpit} estimates.}

\description{

This method computes a matrix of \code{n.boot} subsampled resamples of the regression model. If \code{m} is equal to the dataset size, then a standard bootstrap is performed. Otherwise \code{m} of \code{n} observations, drawn without replacement, are used to re-estimate the model. If \code{m} is null, than a stochastic \code{m} is selected based on a uniform \eqn{\alpha}{alpha} and \eqn{m=n^\alpha}{m=n^alpha}. The minimum value of \eqn{\alpha}{alpha} is chosen so that the minimum m is large enough that the expected number of events at the minimum subsample size is 10 per model parameter.

}

\section{Methods}{

\describe{

\item{subresample}{\code{signature(object = "blm",m=NULL,n.boot=25)}:
Computes a subsample resample of the regression coefficients.
  }

 \item{subresample}{\code{signature(object = "lexpit",m=NULL,n.boot=25)}:
Computes a subsample resample of the regression coefficients.
  }

 }
}



\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\keyword{methods}
