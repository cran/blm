\name{ci.logit.transform}
\docType{methods}
\alias{ci.logit.transform}
\alias{ci.logit.transform,blm-method}
\alias{ci.logit.transform,lexpit-method}

\title{Compute confidence interval for linear combination of estimates from a \code{blm} and \code{lexpit} fit.}

\description{

Computes a confidence interval for the specified linear combination of the regression coefficients of a \code{blm} or \code{lexpit} model using a logit-transform method.

}

\section{Methods}{
\describe{

  \item{ci.logit.transform}{\code{signature(object,C,alpha=.05,sig=4,coef=TRUE,average)}:
Computes an estimate and confidence interval for given linear combination of the model parameters specified by \code{C}, a vector of the values for each covariate of the \code{blm} model. The method is intended for use when the \code{blm} coefficients are actively constrained. A logit-transform is used to construct a confidence region. The argument \code{coef} specifies whether a confidence region for a model parameter is to be obtained or whether a risk is being estimated. The argument \code{average} is a vector of logical values indicating which covariates are to be averaged over (\code{TRUE}) and which are to be treated as fixed at \code{C} (\code{FALSE}).
  }

  \item{ci.logit.transform}{\code{signature(object,C,alpha=.05,sig=4,coef=TRUE,average)}:
Computes an estimate and confidence interval for given linear combination of the model parameters specified by \code{C}, a list with components \code{linear} and \code{expit}, which are each vectors of the values for the linear and expit terms of the \code{lexpit} model. The method is intended for use when the \code{lexpit} coefficients are actively constrained. A logit-transform is used to construct a confidence region. The argument \code{coef} specifies whether a confidence region for a model parameter is to be obtained or whether a risk is being estimated. The argument \code{average} is a list of \code{linear} and \code{expit} vectors, which are logical values indicating which covariates are to be averaged over (\code{TRUE}) and which are to be treated as fixed at their \code{C} values (\code{FALSE}).
  }

 }
}

\seealso{
\code{\link{ci}}
}

\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\examples{

data(birthwt)

# ESTIMATE RISK DIFFERENCE ASSOCIATED WITH MOTHER'S SMOKING STATUS, AGE, WEIGHT

fit <- lexpit(low~smoke,low~I(scale(age))+I(scale(lwt)),data=birthwt)

# RISK OF LOW BIRTH WEIGHT OF CHILD BORN TO MOTHER OF AVERAGE AGE, WEIGHT

ci(fit,C=list(linear=1,expit=c(1,0,0)),coef=FALSE,
	average=list(linear=FALSE,expit=rep(FALSE,3)),
	method="logit")

}


\keyword{methods}
