\name{logLik}
\docType{methods}
\alias{logLik.blm}
\alias{logLik,blm-method}
\alias{logLik,lexpit-method}

\title{Log-likelihood of \code{blm} and \code{lexpit} objects.}

\description{

Method to access the log-likelihood of the fitted \code{blm} or \code{lexpit} model.

}

\section{Methods}{
\describe{

  \item{logLik}{\code{signature(object = "blm",...)}:
Extract log-likelihood. Returns object of \code{logLik} class.
  }

  \item{logLik}{\code{signature(object = "lexpit",...)}:
Extract log-likelihood. Returns object of \code{logLik} class.
  }

 }
}

\details{
	The return object is of the \code{logLik} class. This method is registered with the \link{stats4} package and can therefore be used with applicable methods like \link{AIC} and \link{BIC}.
	Note that when weights are used in the model estimation, the \code{logLik} is a pseduo-log-likelihood.
}

\seealso{\link{logLik.lm}}

\author{Stephanie Kovalchik \email{s.a.kovalchik@gmail.com}}

\examples{

data(ccdata)

fit <- lexpit(y~female, y~packyear, data = ccdata,
       			weight = ccdata$w, strata = ccdata$strata)

logLik(fit)

AIC(fit)
}


\keyword{methods}
