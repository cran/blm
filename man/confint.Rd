\name{confint}
\docType{methods}
\alias{confint}
\alias{confint,blm-method}
\alias{confint,lexpit-method}

\title{Confidence intervals for parameters of \code{blm} and \code{lexpit} objects.}

\description{

Return the confidence intervals for specified parameters and confidence level.

}

\section{Methods}{
\describe{

  \item{confint}{\code{signature(object = "blm", parm, level = 0.95,...)}:
Returns confidence interval (at a given \code{level}) for the specified regression parameters.
  }

  \item{confint}{\code{signature(object = "lexpit", parm, level = 0.95,...)}:
Returns confidence interval (at a given \code{level}) for the specified regression parameters.
  }

 }
}


\author{Stephanie Kovalchik \email{s.a.kovalchik@gmail.com}}

\examples{

data(ccdata)

fit <- lexpit(y~female, y~packyear, data = ccdata,
       			weight = ccdata$w, strata = ccdata$strata)

confint(fit)
}


\keyword{methods}
