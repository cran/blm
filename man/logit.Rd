\name{logit}
\alias{logit}


\title{Logit function}

\description{
Returns the logit. Where,

\deqn{
logit(x) = \log(x/(1-x))
}
}

\usage{
logit(x)
}

\arguments{
  \item{x}{numeric vector}
}

\value{
Numeric that is the logit of \code{x}.
}

\seealso{\code{\link{expit}}}

\examples{

logit(1:10)

}