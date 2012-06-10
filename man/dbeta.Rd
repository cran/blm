\name{leverage}
\docType{methods}
\alias{leverage}
\alias{leverage,blm-method}
\alias{leverage,lexpit-method}

\title{Measures of leverage of observations in a \code{blm} and \code{lexpit} model fit.}

\description{

Returns the leverage of observations in a \code{blm} and \code{lexpit} model fit.

}

\section{Methods}{
\describe{

  \item{leverage}{\code{signature(object = "blm")}:
Extract leverages for the fitted \code{blm} model.

  }

  \item{leverage}{\code{signature(object = "lexpit")}:
Extract leverages for the fitted \code{lexpit} model. 
  }

 }
}


\author{Stephanie Kovalchik \email{s.a.kovalchik@gmail.com}}

\examples{

data(ccdata)

fit <- blm(y~female+packyear, data = ccdata,
       			weight = ccdata$w, strata = ccdata$strata)

summary(leverage(fit))

fit <- lexpit(y~female, y~packyear, data = ccdata,
       			weight = ccdata$w, strata = ccdata$strata)

sort(leverage(fit), decreas=TRUE)[1:5]
}


\keyword{methods}
