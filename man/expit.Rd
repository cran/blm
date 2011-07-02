\name{expit}
\alias{expit}


\title{Inverse-logit function}

\description{
Returns the inverse logit. Where,

\deqn{
expit(x) = \frac{\exp(x)}{(1+\exp(x))}
}
}

\usage{
expit(x)
}

\arguments{
  \item{x}{numeric vector}
}

\value{
Numeric that is the inverse logit of \code{x}.
}

\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\examples{

expit(1:10)

}