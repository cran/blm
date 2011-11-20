\name{ci}
\docType{methods}
\alias{ci}
\alias{ci,blm-method}
\alias{ci,lexpit-method}

\title{Compute confidence interval for linear combination of estimates from a \code{blm} and \code{lexpit} fit.}

\description{

Computes a confidence interval for the specified linear combination of the regression coefficients of a \code{blm} or \code{lexpit} model.

}

\section{Methods}{object,C,alpha=.05,
                   sig=4,
                   method=c("hessian","subsample","logit"),...),
                  
\describe{

  \item{ci}{\code{signature(object = "blm",C="numeric",alpha=0.05,sig=4,method=c("hessian","subsample","logit"),...)}:
Computes an estimate and confidence interval for given linear combination of the model parameters specified by \code{C}. The \code{method} specifies the type of confidence interval to be computed: a large-sample normal (\code{hessian}), a subsampler (\code{subsample}), or a logit-transform based CI.
  }

  \item{ci}{\code{signature(object = "lexpit",C="numeric",alpha=0.05,sig=4,method=c("hessian","subsample","logit"),...)}:
Computes an estimate and confidence interval for given linear combination of the model parameters specified by \code{C}, al list of \code{linear} and \code{expit} components.
  }

Additional arguments \code{...} are passed to \code{ci.blm} or \code{ci.logit.transform} depending on the selected method.
 }
}

\seealso{
\code{\link{subresample}},\code{\link{ci.blm}},\code{\link{ci.logit.transform}}
}

\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\examples{

data(birthwt)

# ADDITIVE MODEL FOR LOW BIRTHWEIGHT
# ESTIMATE RISK DIFFERENCE ASSOCIATED WITH MOTHER'S SMOKING STATUS AND AGE

fit <- blm(low~smoke+I(scale(age)),data=birthwt)

# ESTIMATE RISK DIFFERENCE ASSOCIATED WITH SMOKING

ci(fit,c(0,1,0),method="hessian")

}


\keyword{methods}
