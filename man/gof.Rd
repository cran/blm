\name{gof}
\docType{methods}
\alias{gof}
\alias{gof,blm-method}
\alias{gof,lexpit-method}

\title{Get goodness-of-fit statistic \code{blm} and \code{lexpit} objects.}

\description{

A Hosmer-Lemeshow type goodness-of-fit for \code{blm} or \code{lexpit} model. This test is appropriate when the replicates in each covariate class are few.

}

\section{Methods}{
\describe{

 \item{gof}{\code{signature(object = "blm",groups=10)}:
Hosmer-Lemeshow goodness-of-fit test where groups are based on order of predicted risk.
  }

 \item{gof}{\code{signature(object = "lexpit",groups=10)}:
Hosmer-Lemeshow goodness-of-fit test where groups are based on order of predicted risk.
  }

 }
}

\references{
Hosmer, D.W., Lemeshow, S. (2000) Applied Logistic Regression, New York : Wiley, ISBN 0471615536.
}

\examples{

data(birthwt)

fit.linear <- blm(low~smoke+ht+ui+I(ptl>0),birthwt)

fit.linear

gof(fit.linear)

}

\seealso{\code{\link{predict}},\code{\link{dispersion}}}

\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\keyword{methods}
