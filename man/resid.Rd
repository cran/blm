\name{resid}
\docType{methods}
\alias{resid}
\alias{resid,blm-method}
\alias{resid,lexpit-method}

\title{Get residuals for \code{blm} and \code{lexpit} objects.}

\description{

Computes vector of residuals: the observed event outcome minus the predicted probability.

}

\section{Methods}{
\describe{

 \item{resid}{\code{signature(object = "blm")}:
Residuals for observations of model matrix. 
  }

 \item{resid}{\code{signature(object = "lexpit")}:
Residuals for observations of model matrix.
  }

 }
}


\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\examples{

data(birthwt)

fit.linear <- blm(low~smoke+factor(race)+I(scale(age)),birthwt)

head(resid(fit.linear))

}


\keyword{methods}
