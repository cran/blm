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
    \item{\code{fit}:}{List of output from \code{constrOptim.nl} optimization.}
    \item{\code{par.start}:}{List of \code{par.start} and \code{not.feasible}.}
    \item{\code{f.loglik}:}{Function, BLM log-likelihood (objective function of optimization) that takes single argument that is the vector of model paramters.}
    \item{\code{f.score}:}{Function, BLM score (first derivative of log-likelihood) that takes single argument that is the vector of model paramteres.}
    \item{\code{run.time}:}{Numeric, elapsed time of maximization in seconds.}
    \item{\code{data}:}{Data frame containing the outcomes and explanatory variables.}
    \item{\code{formula}:}{Model formula.}
    \item{\code{ineq}:}{matrix of inequality constraints.}
    \item{\code{n.missing}:}{number of missing observations removed.}
    \item{\code{constraints}:}{List of \code{ineq} and \code{ineq.jac} which specify the set of linear inequalities and their Jacobian. Each are functions that take parameters as their only argument.}
\item{\code{active.constraints}:}{List of final estimates evaluated at \code{ineq} that are less than the \code{tol} setting for \code{blm}. This is the criterion for judging a constraint active. The components of the list are \code{active.constraints}, a matrix of the linear combination of the model parameters at the boundary of the feasible region, \code{boundary.value} the specific values at the value corresponding to the active constraints,
 and \code{active.grad} the Jacobian of active constraints.}
     \item{\code{H}:}{Hessian matrix.}
    \item{\code{V}:}{Taylor series unconstrained variance-covariance matrix.}
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
Extractor for MLEs. 
  }
  
  \item{vcov}{\code{signature(object = "blm")}:
Extractor for variance-covariance based on Taylor series large-sample Hessian approximation.
  }
  
    \item{vcovBoot}{\code{signature(object = "blm",n.boot=25)}:
Computes a bootstrap variance-covariance using \code{n.boot} resamplings with replacement.
  }
  
 \item{gof}{\code{signature(object = "blm",groups=10)}:
Hosmer-Lemeshow goodness-of-fit test where groups are based on ordered rank of predicted risks.
  }

  \item{ci}{\code{signature(C,object = "blm",alpha=0.05,sig=4,...)}:
Computes an estimate and confidence interval for given linear combination of the model parameters specified by code{C}.
  }
 }
}


\keyword{classes}
