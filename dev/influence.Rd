\name{influence}
\docType{methods}
\alias{influence}
\alias{influence,blm-method}
\alias{influence,lexpit-method}

\title{Measures of influence for \code{blm} and \code{lexpit} objects.}

\description{

Returns a matrix \code{influence} and univariate measure \code{summary.influence} which are the score contributions and proportionate squared score, respectively. Larger values of \code{summary.influence} are more suggestive of an influential point.

}

\section{Methods}{
\describe{

 \item{influence}{\code{signature(model = "blm")}:
List of \code{influence} and \code{summary.influence}.
  }

 \item{influence}{\code{signature(model = "lexpit")}:
List of \code{influence} and \code{summary.influence}.
  }

 }
}


\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}


\keyword{methods}
