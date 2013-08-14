\name{LRT}
\alias{LRT}


\title{Performs likelihood-ratio test for lexpit and BLM models of cohort data}


\description{
Computes the likelihood ratio test for the significance of the specified variable in a lexpit or BLM model fit to cohort data. This method is only valid for study designs that use simple random sampling.

}

\usage{
LRT(object, var)
}

\arguments{
  \item{object}{a model of the \code{lexpit} or \code{blm} class.}
  \item{var}{character name of \code{term.label} to be tested}
}

\value{
A matrix with the LRT statistic and p-value for the test of the significance of the specified variable given all other variables in the model.
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

LRT(fit, "x1")
LRT(fit, "x2")

}
