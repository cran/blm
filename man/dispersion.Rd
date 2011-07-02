\name{dispersion}
\docType{methods}
\alias{dispersion}
\alias{dispersion,blm-method}
\alias{dispersion,lexpit-method}

\title{Dispersion statistics for \code{blm} and \code{lexpit} objects.}

\description{

Computes the deviance and Pearson chi-squared statistics for the fit from a \code{blm} or \code{lexpit} model. These tests are appropriate when all predictors are categorical and there are many replicates within each covariate class.

}

\section{Methods}{
\describe{

 \item{dispersion}{\code{signature(object = "blm")}
  }

 \item{dispersion}{\code{signature(object = "lexpit")}
  }
 }
}


\references{
McCullagh, Peter; Nelder, John (1989). Generalized Linear Models, Second Edition. Boca Raton: Chapman and Hall/CRC. ISBN 0-412-31760-5.
}

\examples{

data(birthwt)

fit.linear <- blm(low~smoke+ht+ui+I(ptl>0),birthwt)

fit.linear

dispersion(fit.linear)

}

\seealso{\code{\link{gof}}}

\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\keyword{methods}
