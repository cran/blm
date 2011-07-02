\name{vcovBoot}
\docType{methods}
\alias{vcovBoot}
\alias{vcovBoot,blm-method}
\alias{vcovBoot,lexpit-method}

\title{Get bootstrap variance-covariance from \code{blm} and \code{lexpit} objects.}

\description{

Compute a bootstrap variance-covariance matrix for the coefficients for the fit of a \code{blm} or \code{lexpit} model.

}

\section{Methods}{
\describe{

    \item{vcovBoot}{\code{signature(object = "blm",n.boot=25)}:
Computes a bootstrap variance-covariance using \code{n.boot} resamplings with replacement.
  }
 }
}


\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\keyword{methods}
