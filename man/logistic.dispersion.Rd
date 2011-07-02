\name{logistic.dispersion}
\alias{logistic.dispersion}


\title{Computes deviance and Pearson's chi-squared statistica for a logisitc model fit with \code{glm}.}

\description{

Computes the deviance and Pearson chi-squared statistics for the fit from a \code{glm} logistic regression. These tests are appropriate when all predictors are categorical and there are many replicates within each covariate class.

}

\usage{
logistic.dispersion(f,data,fit)
}

\arguments{
  \item{f}{model formula}
  \item{data}{data frame on which the \code{fit} was based}
  \item{fit}{glm object for logistic regression}
}

\value{
Returns a list of \code{observed}, \code{expected}, \code{pearson}. \code{pearson.p}, \code{pearson.df}, \code{deviance}, \code{deviance.p} and \code{deviance.df}.
}

\examples{

data(birthwt)

fit <- glm(low~smoke+ht+ui+I(ptl>0),birthwt,family=binomial)

fit

logistic.dispersion(low~smoke+ht+ui+I(ptl>0),birthwt,fit)


}

\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}
