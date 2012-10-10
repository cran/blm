\name{gof}
\alias{gof}

\title{Hosmer-lemeshow goodness-of-fit statistics for \code{blm} and \code{lexpit} objects.}

\description{

Computes the deviance and Pearson chi-squared statistics for the fit from a \code{blm} or \code{lexpit} model. These tests are appropriate when all predictors are categorical and there are many replicates within each covariate class.

}

\section{usage}{
gof(object)
}

\section{arguments}{

\describe{

 \item{object}{instance of \code{blm} or \code{lexpit}}
  }
 }

\value{
Returns a list with \code{table}, with expected \code{E} and observed \code{O}, and the chi-square test \code{chisq} and p-value (\code{p.value}) for the Pearson goodness-of-fit test. The observed and expected count are listed in the order of the unique levels formed by the design matrix.

When sample weights are present, the goodness-of-fit test is a modified F-test as suggested by Archer et al. (2007).
}

\examples{

data(ccdata)

ccdata$packyear <- ccdata$packyear+runif(nrow(ccdata))

# UNWEIGHTED GOF
fit <- blm(y~female+packyear,data = ccdata)
gof(fit)

# WEIGHTED GOF
fit <- blm(y~female+packyear,data = ccdata, weight = ccdata$w)
gof(fit)

}

\references{
Archer KJ, Lemeshow S, Hosmer DW. Goodness-of-fit tests for logistic regression models when data are collected using a complex sampling design. \emph{Computational Statistics & Data Analysis}. 2007;51:4450--4464. 
}

\seealso{\code{\link{blm}}, \code{\link{lexpit}}}

\author{Stephanie Kovalchik \email{s.a.kovalchik@gmail.com}}

