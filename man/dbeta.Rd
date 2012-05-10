\name{dBeta}
\docType{methods}
\alias{dBeta}
\alias{dBeta,blm-method}
\alias{dBeta,lexpit-method}

\title{Measures of influence of parameters \code{blm} and \code{lexpit} objects.}

\description{

Returns the jackknife deviations by subject and regression parameter which can be used to assess the influence of rows and columns of the fitted \code{lexpit} or \code{blm} model.

}

\section{Methods}{
\describe{

  \item{dBeta}{\code{signature(object = "blm")}:
Extract dBetas for the fitted \code{blm} model. These are the jackknife deviations by subject and regression coefficient.
  }

  \item{dBeta}{\code{signature(object = "lexpit")}:
Extract dBetas for the fitted \code{lexpit} model. These are the jackknife deviations by subject and regression coefficient.
  }

 }
}


\author{Stephanie Kovalchik \email{s.a.kovalchik@gmail.com}}

\examples{

data(ccdata)

fit <- blm(y~female+packyear, data = ccdata,
       			weight = ccdata$w, strata = ccdata$strata)

apply(dBeta(fit),2,summary) # LOOK AT INFLUENCE RANGE FOR EACH PARAMETER

fit <- lexpit(y~female, y~packyear, data = ccdata,
       			weight = ccdata$w, strata = ccdata$strata)

dBeta(fit)$expit[1:10]
}


\keyword{methods}
