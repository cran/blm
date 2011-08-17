\name{vcov}
\docType{methods}
\alias{vcov}
\alias{vcov,blm-method}
\alias{vcov,lexpit-method}

\title{Get variance-covariance from \code{blm} and \code{lexpit} objects.}

\description{

Returns Hessian-based variance-covariance matrix of the fit of a \code{blm} or \code{lexpit} model. If any constraints are active, only the augmented Lagrangian takes this into account in the Hessian computation, so if \code{augmented} is \code{FALSE}, i.e. the adaptive barrier method of optimization is used, the covariance-variance might be inaccurate.
}

\section{Methods}{
\describe{

\item{vcov}{\code{signature(object = "blm")}:
Extractor for variance-covariance of MLEs. 
  }
 }
}


\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\keyword{methods}
