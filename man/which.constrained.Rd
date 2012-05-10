\name{which.at.boundary}
\alias{which.at.boundary}

\title{Covariate patterns at the boundary for \code{blm} and \code{lexpit} objects.}

\description{

Returns matrix of covariate types with a predicted probability at the lower or upper boundary defined by the specified \code{criterion} or \code{NA} if no boundary constraints are present.

}

\section{usage}{
which.at.boundary  (object, criterion = 1e-06) 
}

\section{arguments}{
\describe{
\item{object}{model fit of class \code{blm} or \code{lexpit}}
\item{criterion}{numeric distance from 0 (or 1) that is considered to be at the boundary}
 }
}

\value{
Returns all rows of design matrix whose predicted risk are less than or equal to \code{criterion} or greater than or equal to 1 - \code{criterion}.
}

\author{Stephanie Kovalchik \email{s.a.kovalchik@gmail.com}}

\examples{

data(ccdata)

fit <- blm(y~female+packyear,data = ccdata,
			weight = ccdata$w, strata = ccdata$strata)

which.at.boundary(fit)

}
