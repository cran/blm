\name{blm-class}
\docType{class}
\alias{blm-class}

\title{Class "blm"}

\description{
Class for binomial linear regression (BLM).
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("blm", ...)}. 
}

\section{Slots}{
  \describe{
  \item{\code{coef}:}{vector of fitted coefficients}
  \item{\code{vcov}:}{matrix of variance-covariate estimates for \code{coef}}
  \item{\code{formula}:}{model formula}
  \item{\code{df.residual}:}{residual degrees of freedom}
  \item{\code{data}:}{data frame used in fitting, after applying \code{na.action}}
  \item{\code{which.kept}:}{vector of index of values in original data source that were used in the model fitting}
  \item{\code{y}:}{response vector for fitted model}	
  \item{\code{weights}:}{vector of weights used in model fitting}
  \item{\code{strata}:}{stratification factor for weighted regression.}
  \item{\code{converged}:}{logical message about convergence status at the end of algorithm}
  \item{\code{par.init}:}{initial parameter values for optimization algorithm}
  \item{\code{loglik}}{value of log-likelihood (normalized for weighted likelihood) under full model} 
  \item{\code{loglik.null}}{value of log-likelihood (normalized for weighted likelihood) under null model} 
  \item{\code{barrier.value}}{value of the barrier function at the optimum}
   \item{\code{dBeta}}{matrix of influence measures}
 } 
}			     


\section{Methods}{

  \describe{

\item{show}{\code{signature(object = "blm")}:
Display point estimates of \code{blm} object.
  }

\item{print}{\code{signature(x = "blm",...)}:
Display point estimates of \code{blm} object.
  }

\item{summary}{\code{signature(object = "blm",...)}:
List of estimates and convergence information.
  }

\item{coef}{\code{signature(object = "blm")}:
Extractor for fitted coefficients. 
  }

\item{resid}{\code{signature(object = "blm")}:
Extractor for residuals. 
  }
  
  
  \item{vcov}{\code{signature(object = "blm")}:
Extractor for variance-covariance based on Taylor series large-sample Hessian approximation with the pseudo-likelihood of the constrained optimization.
  }

  \item{predict}{\code{signature(object = "blm")}:
Returns vector of linear predictors for each subject of the fitted model.
  }

  \item{confint}{\code{signature(object = "blm", parm, level = 0.95,...)}:
Returns confidence interval (at a given \code{level}) for the specified regression parameters.
  }

  \item{dBeta}{\code{signature(object = "blm")}:
Extract dBetas for the fitted \code{blm} model. These are the jackknife deviations by subject and regression coefficient.
  }

 }
}

\seealso{\code{\link{blm}}, \code{\link{constrOptim}}}

\keyword{classes}
