\name{EO}
\docType{methods}
\alias{EO}
\alias{EO,blm-method}
\alias{EO,lexpit-method}

\title{Compute the ratio of expected event to observed events for \code{blm} and \code{lexpit} objects.}

\description{

Returns a list of expected to observed counts and the specified confidence interval. The argument \code{group} can be used to estimate this ratio by the categories of the categorical variable \code{group}.

}

\section{Methods}{
\describe{

 \item{EO}{\code{signature(object = "blm",alpha=.05,group=NULL,...)}:
Expected to observed based on probability predictions from the \code{blm} fit.
  }

 \item{EO}{\code{signature(object = "lexpit",alpha=.05,group=NULL,...)}:
Expected to observed based on probability predictions from the \code{lexpit} fit.

  }

 }
}


\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\examples{

data(birthwt)

fit.linear <- blm(low~smoke+factor(race)+I(scale(age)),birthwt)

EO(fit.linear)

}


\keyword{methods}
