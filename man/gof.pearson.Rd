\name{gof.pearson}
\alias{gof.pearson}

\title{Pearson's goodness-of-fit statistics for \code{blm} and \code{lexpit} objects.}

\description{

Computes the deviance and Pearson chi-squared statistics for the fit from a \code{blm} or \code{lexpit} model. These tests are appropriate when all predictors are categorical and there are many replicates within each covariate class.

}

\section{usage}{
gof.pearson(object)
}

\section{arguments}{

\describe{

 \item{object}{instance of \code{blm} or \code{lexpit}}
  }
 }

\value{
Returns a list with expected \code{E} and observed \code{O} and the chi-square test \code{chisq} and p-value (\code{p.value}) for the Pearson goodness-of-fit test. The observed and expected count are listed in the order of the unique levels formed by the design matrix.
}

\examples{

data(ccdata)

fit <- blm(y~female+I(packyear>20),data = ccdata,
			weight = ccdata$w, strata = ccdata$strata)

gof.pearson(fit)

}


\seealso{\code{\link{blm}}, \code{\link{lexpit}}}

\author{Stephanie Kovalchik \email{s.a.kovalchik@gmail.com}}

