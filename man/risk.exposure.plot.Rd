\name{risk.exposure.plot}
\alias{risk.exposure.plot}


\title{Risk-exposure scatter plot}

\description{

Calculates the weighted average crude risk against the average exposure level for a continuous exposure. Each point corresponds to overlapping subgroups of 20 percent of the sample ordered from lowest to highest exposure and a sliding window of 1% of the sample.

}

\usage{
risk.exposure.plot(object, scale=1,...)
}

\arguments{
  \item{object}{list or data.frame with \code{risk} and \code{x} covariate. Return object of \link{crude.risk}}
  \item{scale}{multiplicative factor to modify scale of crude risk estimates}
  \item{...}{additional arguments passed to \link{scatter.smooth}}
}

\details{
The risk-exposure scatter plot is intended to explore the possible functional relationship between risk and exposure.

}

\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\examples{

data(aarp)

risk <- crude.risk(bladder70~redmeat,
		  weights = aarp$w,
		  data = aarp)

risk.exposure.plot(risk,
		   xlab = "Avg. Red Meat Consumption")

}
