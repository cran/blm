\name{plot}
\docType{methods}
\alias{plot}
\alias{plot,blm-method}
\alias{plot,lexpit-method}

\title{Plot of observed versus expected for  \code{blm} and \code{lexpit} objects.}

\description{

Plot the expected versus observed count for a model with fewer than 10 covariate categories, or a predictiveness plot when there are 10 or more categories. The predictiveness plot gives the risk prediction on the y-axis, the risk percentile in the dataset on the x-axis, a red line denoting the observed mean event rate in the sample, and points indicated the observed counts at the midpoint of the decile risk groups (equivalent to the observed groups in a Hosmer-Lemeshow goodness-of-fit test).

}

\section{Methods}{
\describe{

 \item{plot}{\code{signature(x = "blm")}:
Goodness-of-fit plot.
  }

 \item{plot}{\code{signature(x = "lexpit")}:
Goodness-of-fit plot.
  }

 }
}


\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}


\keyword{methods}
