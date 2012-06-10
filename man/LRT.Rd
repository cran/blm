\name{LRT}
\alias{LRT}


\title{Performs likelihood-ratio test for LEXPIT and BLM models of cohort data}


\description{


}

\usage{
LRT(object, ...)
}

\arguments{
  \item{object}{a model of the \code{lexpit} or \code{blm} class.}
  \item{...}{Additional arguments passed to \code{constrOptim}}
}

\value{
A matrix with the LRT statistic and p-value for the test of the significance of each value given all other variables in the model.
}


\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\seealso{\link{constrOptim}}

\examples{


cohort <- data.frame(
	x1 = runif(500),
	x2 = runif(500)
)

cohort$event <- rbinom(n=nrow(cohort),size=1,
			prob=0.25+0.1*cohort$x1+.1*cohort$x2)

fit <- blm(event~x1+x2, data=cohort)

summary(fit)

LRT(fit)

}