\name{lexpit-class}
\docType{class}
\alias{lexpit-class}

\title{Class "lexpit"}

\description{
Class for linear-expit regression (lexpit).
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("lexpit", ...)}. 
}

\section{Slots}{
  \describe{
  \item{\code{coef.linear}:}{vector of fitted linear coefficients}
  \item{\code{coef.expit}:}{vector of fitted expit coefficients}
  \item{\code{vcov.linear}:}{matrix of variance-covariate estimates for linear \code{coef}}
  \item{\code{vcov.expit}:}{matrix of variance-covariate estimates for expit \code{coef}}
  \item{\code{formula.linear}:}{model formula for linear component}
  \item{\code{formula.expit}:}{model formula for expit component}
  \item{\code{df.residual}:}{residual degrees of freedom}
  \item{\code{p}:}{number of linear parameters}
  \item{\code{q}:}{number of expit parameters}
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
   \item{\code{control.lexpit}}{list with control parameters for optimization algorithm}
 } 
}			     


\section{Methods}{

  \describe{

\item{show}{\code{signature(object = "lexpit")}:
Display point estimates of \code{lexpit} object.
  }

\item{print}{\code{signature(x = "lexpit",...)}:
Display point estimates of \code{lexpit} object.
  }

\item{summary}{\code{signature(object = "lexpit",...)}:
List of estimates and convergence information.
  }

\item{coef}{\code{signature(object = "lexpit")}:
Extractor for fitted coefficients. 
  }

\item{logLik}{\code{signature(object = "lexpit")}:
Extractor for log-likelihood of \code{lexpit} model.
  }

\item{model.formula}{\code{signature(object = "lexpit")}:
Extractor for formula of \code{lexpit} object. 
  }

  \item{vcov}{\code{signature(object = "lexpit")}:
Extractor for variance-covariance based on Taylor series large-sample Hessian approximation with the pseudo-likelihood of the constrained optimization.
  }

\item{resid}{\code{signature(object = "lexpit")}:
Extractor for residuals. 
  }

  \item{predict}{\code{signature(object = "lexpit")}:
Returns vector of linear predictors for each subject of the fitted model.
  }

  \item{confint}{\code{signature(object = "lexpit", parm, level = 0.95,...)}:
Returns confidence interval (at a given \code{level}) for the specified regression parameters.
  }

 }
}

\seealso{\code{\link{lexpit}}, \code{\link{constrOptim}}}

\keyword{classes}
