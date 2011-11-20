\name{which.constrained}
\docType{methods}
\alias{which.constrained}
\alias{which.constrained,blm-method}
\alias{which.constrained,lexpit-method}

\title{Index of constrained observations for \code{blm} and \code{lexpit} objects.}

\description{

Returns list of the row index for constrained observations and a vector of the total influence of the constrained covariate profile. The \code{influence}/m can be interpreted as the excess representation of the observation in the estimating equations, where m is the total number of observations with the same covariate pattern.

}

\section{Methods}{
\describe{

 \item{which.constrained}{\code{signature(object = "blm")}:
List of index and influence of constrained observations.
  }

 \item{which.constrained}{\code{signature(object = "lexpit")}:
List of index and influence of constrained observations.
  }

 }
}


\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}


\keyword{methods}
