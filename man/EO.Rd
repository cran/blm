\name{EO}
\alias{EO}

\title{Compute the ratio of expected event to observed events for \code{blm} and \code{lexpit} objects.}

\description{

Returns a list of expected to observed counts and the specified confidence interval. The argument \code{group} can be used to estimate this ratio by the categories of the categorical variable \code{group}. If population-based case-control data is used to fit the model, the expected counts are for the population and make use of the sampling weights. 

}

\usage{
EO(object, index = NULL, level = 0.95)
}

\arguments{
\item{object}{object of class \code{blm} or \code{lexpit}}
\item{index}{factor for computing E/O comparison by subgroups}
\item{level}{numeric, confidence level (between 0 and 1) for the E/O ratios}
}

\value{
Data frame with:
\describe{

\item{\code{E}}{expected count}
\item{\code{O}}{observed counts}
\item{\code{EtoO}}{ratio of expected to observed}
\item{\code{lowerCI}}{lower endpoint of confidence interval for E over O ratio}
\item{\code{upperCI}}{upper endpoint of confidence interval for E over O ratio}
}
}
\author{Stephanie Kovalchik \email{s.a.kovalchik@gmail.com}}

\examples{

data(ccdata)

fit <- blm(y~female+packyear,data = ccdata,
			weight = ccdata$w, strata = ccdata$strata)

EO(fit)

EO(fit, ccdata$strata) # BY FACTOR
}

