\name{ci.blm}
\docType{methods}
\alias{ci.blm}
\alias{ci.blm,blm-method}
\alias{ci.blm,lexpit-method}

\title{Compute confidence interval for linear combination of estimates from a \code{blm} and \code{lexpit} fit.}

\description{

Computes a confidence interval for the specified linear combination of the regression coefficients of a \code{blm} or \code{lexpit} model.

}

\section{Methods}{
\describe{

  \item{ci.blm}{\code{signature(object = "blm",C="numeric",alpha=0.05,sig=4,subsample=FALSE,...)}:
Computes an estimate and confidence interval for given linear combination of the model parameters specified by \code{C}. The \code{subsample} argument requests confidence intervals based on a resampled subsample with \code{...} providing additional arguments to \code{subresample}.
  }

  \item{ci.blm}{\code{signature(object = "lexpit",C="list",alpha=0.05,sig=4,baseline=TRUE,subsample=FALSE,...)}:
Computes an estimate and confidence interval for given linear combination of the model parameters specified by \code{C}, a list with named components \code{linear} and \code{expit}. If \code{baseline} is included, then the list argument \code{expit} supplies the values for the expit components, including the intercept. The \code{subsample} argument requests confidence intervals based on a resampled subsample with \code{...} providing additional arguments to \code{subresample}.
  }

 }
}

\seealso{
\code{\link{subresample}},\code{\link{ci}}
}

\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\examples{

data(birthwt)

# ADDITIVE MODEL FOR LOW BIRTHWEIGHT
# ESTIMATE RISK DIFFERENCE ASSOCIATED WITH MOTHER'S SMOKING STATUS AND AGE

fit <- blm(low~smoke+I(scale(age)),data=birthwt)

# ESTIMATE RISK DIFFERENCE ASSOCIATED WITH SMOKING

ci.blm(fit,c(0,1,0))

}


\keyword{methods}
