\name{aarp}
\docType{data}
\alias{aarp}

\title{Nested case-control data set of bladder cancer in the NIH-AARP Diet and Health Study}

\description{The \code{aarp} data set is a nested case-control study of bladder cancer outcomes in the NIH-AARP Diet and Health Study. The data set is intended for demonstration purposes only.}

\usage{
aarp
}
\format{
  \tabular{ll}{
    bladder70:\tab indicator of bladder cancer by age 70 years \cr
    female:\tab indicator of female gender \cr
    smoke_status: \tab factor of smoking status (four categories)  \cr
    w: \tab inverse of sampling fraction\cr
    redmeat: \tab total daily redmeat consumption (grams/day)\cr
    fiber.centered: \tab total daily fiber consumption (grams), centered on sample median \cr
    educ: \tab factor of education status (six categories)
  }
}

\source{
National Cancer Institute. National Institutes of Health AARP Diet and Health Study. \url{https://dceg.cancer.gov/research/who-we-study/nih-aarp-diet-health-study}. Accessed: 12/10/2012
}

\examples{

data(aarp)

# ABSOLUTE RISK OF BLADDER CANCER BY 70 YEARS
# FOR DIFFERENT GENDER AND RISK GROUP

fit <- blm(bladder70~female * smoke_status, 
			      data = aarp, 
			      weight=aarp$w)

# INTERCEPT IS BASELINE RISK
# ALL OTHER COEFFICIENTS ARE RISK DIFFERENCES FROM BASELINE

summary(fit)


}

\keyword{datasets}
