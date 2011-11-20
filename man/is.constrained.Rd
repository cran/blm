\name{is.constrained}
\docType{methods}
\alias{is.constrained}
\alias{is.constrained,blm-method}
\alias{is.constrained,lexpit-method}

\title{Indicator whether final estimate of \code{blm} or \code{lexpit} fit was constrained.}

\description{

Returns \code{TRUE} or \code{FALSE}. If \code{TRUE}, at least one estimated risk is at the boundary.

}

\section{Methods}{
\describe{

 \item{is.constrained}{\code{signature(object = "blm")}:
Logical
  }

 \item{is.constrained}{\code{signature(object = "lexpit")}:
Logical
  }

 }
}


\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}


\keyword{methods}
