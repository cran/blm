\name{excess.risk}
\alias{excess.risk}


\title{Compute binned excess risk for BLM of LEXPIT fit}

\description{

Calculates the weighted average of the residual (excess) risk for a BLM or LEXPIT model by a specified grouping variable.

}

\usage{
excess.risk(object, group, linear.only=FALSE)
}

\arguments{
  \item{object}{fitted \code{\link{blm}} of \code{\link{lexpit}} model}
  \item{group}{vector, grouping variable to bin excess risk on. Should be the same length as the design matrix for the model fit (after applying na.action).}
  \item{linear.only}{logical, indicator for LEXPIT objects whether excess risk should only adjust for linear component of model.}
}

\value{
Vector of weighted average of residual risks whose names are the factor groups.
}

\details{
The excess risk can be used to assess the functional relationship of exposures not included in the fitted model. 

}

\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\examples{

data(ccdata)

# INVESTIGATE THE FUNCTIONAL RELATIONSHIP BETWEEN RISK AND PACKYEARS
# AFTER ACCOUNTING FOR GENDER

fit <- blm(y~female, weights = ccdata$w,
       			     	strata=ccdata$strata,
				data=ccdata)

r <- excess.risk(fit, group = ccdata$packyear)
r

plot(y = r, x = as.numeric(names(r)),
       ylab = "Binned excess risk",
       xlab = "packyears", las = 1,
       type = "b")

}
