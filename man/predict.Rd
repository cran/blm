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


\author{Stephanie Kovalchik \email{s.a.kovalchik@gmail.com}}

\examples{

data(ccdata)

fit <- lexpit(y~female, y~packyear, data = ccdata,
       			weight = ccdata$w, strata = ccdata$strata)

predict(fit)[1:10]
}


\keyword{methods}
