\name{displacement}
\docType{methods}
\alias{displacement}
\alias{displacement,blm-method}
\alias{displacement,lexpit-method}

\title{Measures of confidence interval displacement of observations in a \code{blm} and \code{lexpit} model fit.}

\description{

Returns the confidence interval displacement of observations in a \code{blm} and \code{lexpit} model fit, which measures the influence of each observation on the regression estimates. Large values are suggestive of an influential point.

}

\section{Methods}{
\describe{

  \item{displacement}{\code{signature(object = "blm")}:
Extract displacements for the fitted \code{blm} model.

  }

  \item{displacement}{\code{signature(object = "lexpit")}:
Extract displacements for the fitted \code{lexpit} model. 
  }

 }
}


\author{Stephanie Kovalchik \email{s.a.kovalchik@gmail.com}}

\examples{

cohort <- data.frame(
	x1 = runif(100),
	x2 = runif(100)
)

cohort$event <- rbinom(n=nrow(cohort),size=1,
			prob=0.25+0.1*cohort$x1+.1*cohort$x2)

fit <- blm(event~x1+x2, data=cohort)

C <- displacement(fit)

# INDEX PLOT
plot(C)

}


\keyword{methods}
