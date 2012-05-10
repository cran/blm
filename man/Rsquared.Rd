\name{Rsquared}
\alias{Rsquared}

\title{Compute R-squared measures of model fit for \code{blm} and \code{lexpit} objects.}

\description{

Returns McFadden's unadjusted and adjusted R-squared measures for models of a binary outcome. 
}

\usage{
Rsquared(object)
}

\arguments{
\item{object}{object of class \code{blm} or \code{lexpit}}
}

\value{
List of \code{R2} and \code{R2adj}.
}

\author{Stephanie Kovalchik \email{s.a.kovalchik@gmail.com}}

\examples{

data(ccdata)

fit <- blm(y~female+packyear,data = ccdata,
			weight = ccdata$w, strata = ccdata$strata)

Rsquared(fit)
}

