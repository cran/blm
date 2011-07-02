\name{print}
\docType{methods}
\alias{print}
\alias{print,blm-method}
\alias{print,lexpit-method}

\title{Print coefficients of \code{blm} and \code{lexpit} model fit.}

\description{

Prints the regression coefficients of the fit of a \code{blm} or \code{lexpit} model.

}

\section{Methods}{
\describe{
\item{print}{\code{signature(x = "blm",digits=4,...)}:
Display call, summary matrix of estimates, t-values, p-values  (with \code{digits} number of significant digits), degrees of freedom, Log-lik, and AIC of \code{blm} object.
  }

\item{print}{\code{signature(x = "lexpit",digits=4,...)}:
Display call, summary matrix of estimates, t-values, p-values (with \code{digits} number of significant digits), degrees of freedom, Log-lik, and AIC of \code{lexpit} object.
  }
  }
}



\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\keyword{methods}
