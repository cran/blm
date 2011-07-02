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
\item{gradient}{partial derivative (score) evaluated at final estimates}
\item{feasible}{logical, whether all risks predicted from study design matrix are probabilities}
\item{convergence}{scalar code on convergence result (0 indicates successful convergence)}
\item{message}{character message when convergence is not achieved}
\item{loglik}{numeric log-likelihood at final model estimates}
\item{df}{degrees of freedom}
\item{null.deviance}{deviance, twice the negative of an intercept-only fit}
 }
}


\seealso{For further details on convergence report \code{\link{constrOptim.nl}}}

\examples{

data(grad)

admission <- lexpit(admit~I(rank>1),admit~I(scale(gre)),grad)

summary(admission)

}


\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\keyword{methods}
