\name{ci}
\docType{methods}
\alias{ci}
\alias{ci,blm-method}
\alias{ci,lexpit-method}

\title{Compute confidence interval for linear combination of estimates from a \code{blm} and \code{lexpit} fit.}

\description{

Computes a confidence interval for the specified linear combination of the regression coefficients of a \code{blm} or \code{lexpit} model.

}

\section{Methods}{
\describe{

  \item{ci}{\code{signature(object = "blm",C="numeric",alpha=0.05,sig=4,...)}:
Computes an estimate and confidence interval for given linear combination of the model parameters specified by \code{C}.
  }

  \item{ci}{\code{signature(object = "lexpit",C="numeric",alpha=0.05,sig=4,...)}:
Computes an estimate and confidence interval for given linear combination of the model parameters specified by \code{C}.
  }

 }
}


\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\examples{

data(birthwt)

# ADDITIVE MODEL FOR LOW BIRTHWEIGHT
# ESTIMATE RISK DIFFERENCE ASSOCIATED WITH MOTHER'S SMOKING STATUS AND AGE

fit <- blm(low~smoke+I(scale(age)),data=birthwt)

# ESTIMATE RISK DIFFERENCE ASSOCIATED WITH SMOKING

ci(fit,c(0,1,0))

}


\keyword{methods}
