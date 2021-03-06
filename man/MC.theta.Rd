% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StatCompR.R
\name{MC.theta}
\alias{MC.theta}
\title{A function using a Monte Carlo simulation by the antithetic variate approach and by the simple Monte Carlo method}
\usage{
MC.theta(R = 10000, antithetic = TRUE)
}
\arguments{
\item{R}{the number of samples}

\item{antithetic}{the choice of the method}
}
\value{
the estimation
}
\description{
A function using a Monte Carlo simulation by the antithetic variate approach and by the simple Monte Carlo method
}
\examples{
\dontrun{
set.seed(12345)
m=1000
MC1=MC2=numeric(m)
for (i in 1:m) {
 MC1[i] <- MC.theta(R = 1000, anti = FALSE)
 MC2[i] <- MC.theta(R = 1000)
}
}
}
