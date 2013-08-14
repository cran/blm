\name{crude.risk}
\alias{crude.risk}


\title{Risk-exposure scatter plot}

\description{

Calculates the weighted average crude risk against the average exposure level for a continuous exposure. Each point corresponds to overlapping subgroups of 20 percent of the sample ordered from lowest to highest exposure and a sliding window of 1% of the sample.

}

\usage{
crude.risk(formula, data, weights = NULL, na.action = na.omit) 
}

\arguments{
  \item{formula}{formula specifying the binary outcome and the continuous covariate of interest, e.g. \code{y~x}}
  \item{data}{dataframe containing the variables specified in \code{formula}}
  \item{weights}{vector of sample weights}
  \item{na.action}{function used for handling missing variables in the variables of \code{formula} and \code{weights}}
}

\details{
The \code{crude.risk} function is intended to explore the possible functional relationship between risk and exposure in a non-parametric way.

}


\seealso{\link{risk.exposure.plot}}

\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\examples{

data(aarp)

risk <- crude.risk(bladder70~redmeat,
		  weights = aarp$w,
		  data = aarp)

risk.exposure.plot(risk,
		   xlab = "Avg. Red Meat Consumption")

}
