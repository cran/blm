\name{summary}
\docType{methods}
\alias{summary}
\alias{summary,blm-method}
\alias{summary,lexpit-method}

\title{Summary of \code{blm} and \code{lexpit} model fit.}

\description{

A list of estimates and convergence status of a \code{blm} or \code{lexpit} model fit.

}

\section{Methods}{
\describe{
\item{summary}{\code{signature(x = "blm",...)}:
List of estimates and convergence information.
  }

\item{summary}{\code{signature(x = "lexpit",...)}:
List of estimates and convergence information.
  }

The list returned has the named components:

\item{est}{vector of estimated regression coefficients. For lexpit model estimates are split into \code{est.linear} and \code{est.expit} components of list}
\item{gradient}{partial derivative (score) evaluated at final estimates.}
\item{feasible}{logical, whether all risks predicted from study design matrix are probabilities}
\item{active}{matrix, covariate patterns for risk estimates at the boundary of the feasible region. If \code{NULL}, there are no active constraints.}
\item{convergence}{scalar code on convergence result (0 indicates successful convergence)}
\item{message}{character message when convergence is not achieved}
\item{loglik}{numeric log-likelihood at final model estimates}
\item{df}{degrees of freedom}
\item{null.deviance}{deviance, twice the negative of an intercept-only fit}
\item{R2.efron}{Efron's (1978) pseudo R-squared}
\item{R2.mcfadden}{McFadden's (1974) pseudo R-squared}
\item{R2.mcfadden.adj}{McFadden's adjusted pseudo R-squared}
\item{R2.coxsnell}{Cox and Snell's (1989) pseudo R-squared}
\item{R2.coxsnell.adj}{Cox and Snell's adjusted R-squared}
\item{seconds.to.run}{user elapsed seconds to fit model}
 }
}


\details{
In unconstrained optimization, the elements of the gradient should be near to zero at the maximizing value of the parameters. However, when constraints are active, i.e. when the parameters are at the boundary of the feasible region, this condition might not be met.

An inequality constraint is declared active if, at the final estimates of the algorithm, the evaluated inequality is near to zero. The default threshold is 1e-4. 
}

\seealso{For further details on convergence report \code{\link{constrOptim.nl}}}

\references{

Cox, D. R. & Snell, E. J. (1989). \emph{The Analysis of Binary Data}, Second Edition, London: Chapman and Hall.

Efron, B. (1978), Regression and ANOVA with zero-one data: measures of residual variation, \emph{Journal of the American Statistical Association} 73, 113-121.

McFadden, D. (1974), Conditional logit analysis of qualitative choice behaviour, in: P. Zarembka (ed.), \emph{Frontiers in Econometrics}, Academic Press, New York, 105-142.
}

\examples{

data(grad)

admission <- lexpit(admit~I(rank>1),admit~I(scale(gre)),grad)

summary(admission)

}


\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\keyword{methods}
