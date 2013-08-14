\name{model.formula}
\docType{methods}
\alias{model.formula}
\alias{model.formula,blm-method}
\alias{model.formula,lexpit-method}

\title{Get formula call for \code{blm} and \code{lexpit} objects.}

\description{

Extract vector of formula of the fit of a \code{blm} or the formulas for the \code{lexpit} model.

}

\section{Methods}{
\describe{
\item{model.formula}{\code{signature(object = "blm")}:
Extractor for formula of \code{blm} object.
  }
\item{model.formula}{\code{signature(object = "lexpit")}:
Extractor for formulas of \code{lexpit} object. Returns a list containing the \code{linear} and \code{expit} formulas.
  }
 }
}


\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\keyword{methods}
