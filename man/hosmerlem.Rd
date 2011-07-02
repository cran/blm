\name{hosmerlem}
\alias{hosmerlem}


\title{Hosmer-Lemeshow goodness-of-fit for logistic regression model}

\description{
Computes the chi-squared statistic based on \code{groups} for a logistic regression fit from \code{glm}. The observations are gouped according to the ranked order of the project risks.
}

\usage{
hosmerlem(f,fit,data,groups=10)
}

\arguments{
  \item{f}{formula given to \code{glm}}
   \item{fit}{\code{glm} object from logistic regression fit}
\item{data}{data frame containing the variables of \code{f}}
   \item{groups}{scalar number of groups to bin for the observed - expected deviations}
}

\value{
List of \code{chisq} and \code{p.value}.
}

\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\seealso{\code{\link{glm}}}
