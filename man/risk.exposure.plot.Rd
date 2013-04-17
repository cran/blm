\name{risk.exposure.plot}
\alias{risk.exposure.plot}


\title{Risk-exposure scatter plot}

\description{

Calculates the weighted average crude risk against the average exposure level for a continuous exposure. Each point corresponds to overlapping subgroups of 20 percent of the sample ordered from lowest to highest exposure and a sliding window of 1% of the sample.

}

\usage{
risk.exposure.plot(y, group, weights=NULL,...)
}

\arguments{
  \item{y}{vector of binary event indicators}
  \item{group}{vector of the continuous variable from which the sequential groups will be constructed}
  \item{weights}{vector of sample weights for computing the population-based crude risk}
  \item{...}{additional arguments passed to \link{scatter.smooth}}
}

\value{
Returns a list of \code{risk} and \code{x} containing the average crude risk estimates and covariates, respectively. 
}

\details{
The risk-exposure scatter plot is intended to explore the possible functional relationship between risk and exposure.

}

\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\examples{

data(aarp)

risk.exposure.plot(aarp$bladder70,
		   aarp$redmeat,
		   weights = aarp$w,
		   xlab = "Avg. Red Meat Consumption")

}
