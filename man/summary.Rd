\name{summary}
\docType{methods}
\alias{summary.blm}
\alias{summary,blm-method}
\alias{summary,lexpit-method}

\title{Summary of \code{blm} and \code{lexpit} model fit.}

\description{

A list of estimates and convergence status of a \code{blm} or \code{lexpit} model fit.

}

\section{Methods}{
\describe{
\item{summary}{\code{signature(object = "blm")}:
Matrix of estimates and convergence information.
  }

\item{summary}{\code{signature(object = "lexpit")}:
Matrix of estimates and convergence information.
  }

The matrix returned has the named components:

\item{Est.}{vector of estimated regression coefficients. For lexpit model estimates are split into \code{est.linear} and \code{est.expit} components of list}
\item{Std. Err}{standard error of model estimates}
\item{t-value}{t-value of model estimates}
\item{p-value}{p-value (two-sided) of model estimates}
}
}


\seealso{\code{\link{blm}}, \code{\link{lexpit}}}


\examples{

data(ccdata)

fit <- blm(y~female+packyear,data = ccdata,
			weight = ccdata$w, strata = ccdata$strata)

summary(fit)

fit.lexpit <- lexpit(y~female, y~packyear,data = ccdata,
			weight = ccdata$w, strata = ccdata$strata)

summary(fit.lexpit)


}


\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\keyword{methods}
