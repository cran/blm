\name{show}
\docType{methods}
\alias{show}
\alias{show,blm-method}
\alias{show,lexpit-method}

\title{Show \code{blm} and \code{lexpit} model fit.}

\description{

Print estimates of a \code{blm} or \code{lexpit} model fit.

}

\section{Methods}{
\describe{
\item{show}{\code{signature(x = "blm")}:
Display call, summary matrix of estimates, t-values, p-values  (to 3 significant digits), degrees of freedom, Log-lik, and AIC of \code{blm} object.
  }

\item{show}{\code{signature(x = "lexpit")}:
Display call, summary matrix of estimates, t-values, p-values  (with \code{digits} number of significant digits), degrees of freedom, Log-lik, and AIC of \code{lexpit} object.
  }
  }
 }



\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\keyword{methods}
