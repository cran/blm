\name{ssc}
\alias{ssc}


\title{Strong Sufficient Condition check for binomial linear regression models}

\description{
Check whether SSC are satisfied when \code{blm} or \code{lexpit} have active constraints.
}

\usage{
ssc(object)
}

\arguments{
  \item{object}{object of class \code{blm} or \code{lexpit}}
}

\value{
Returns list of SSC criteria. If each of the following are met then SSC is satisfied and the constrained optimizer converges to a local constrained maximum.

\describe{
\item{barrier.value}{Final value of barrier parameter, should be positive}
\item{row.rank}{Rank of the Jacobian for active constraints. Should be equal to the number of rows}
\item{is.positive.definite}{Check if the matrix product of the Null space of the Jacobian for active constraints is positive definite. The matrix is positive definite if all eigenvalues are positive.}
 }
}

\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}
