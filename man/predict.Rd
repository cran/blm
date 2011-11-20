\name{predict}
\docType{methods}
\alias{predict}
\alias{predict,blm-method}
\alias{predict,lexpit-method}

\title{Get risk predictions for \code{blm} and \code{lexpit} objects.}

\description{

Computes vector of risk predictions for the dataset used to fit the model.

}

\section{Methods}{
\describe{

 \item{predict}{\code{signature(object = "blm")}:
Risk predictions for fit design matrix. 
  }

 \item{predict}{\code{signature(object = "lexpit")}:
Risk predictions for fit design matrix. 
  }

 }
}


\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\examples{

data(birthwt)

fit.linear <- blm(low~smoke+factor(race)+I(scale(age)),birthwt)

head(predict(fit.linear))

}


\keyword{methods}
