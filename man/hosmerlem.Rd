\name{hosmerlem}
\alias{hosmerlem}


\title{Hosmer-Lemeshow goodness-of-fit for binary data model}

\description{
Computes the chi-squared statistic based on \code{groups} for a regression analysis on a binary outcome. The observations are gouped according to the ranked order of the projected risks.
}

\usage{
hosmerlem(y,predicted,groups=10,plot=FALSE,risk.labels=TRUE,sig=3,...)
}

\arguments{
   \item{y}{observed cases}
   \item{predicted}{predicted values}
   \item{groups}{scalar number of groups to bin for the observed - expected deviations}
   \item{plot}{logical whether plot of observed and expected counts should be returned}
   \item{risk.labels}{logical whether labels of midpoint of risk deciles should be included}
   \item{sig}{scalar value for rounding risk labels}
   \item{...}{further arguments to \code{text} to handle rendering of risk labels}
}

\value{
List of \code{chisq}, \code{p.value}, observed counts \code{O}, expected counts \code{E}, table \code{group} showing the quantile groupings.
}

\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\seealso{\code{\link{glm}}}

\examples{

data(grad)

fit <- blm(admit~gpa+gre,data=grad)

hosmerlem(grad$admit,predict(fit))

}