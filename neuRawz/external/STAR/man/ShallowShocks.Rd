\name{ShallowShocks}
\alias{ShallowShocks}
\docType{data}
\title{Shallow Shocks (M >= 6.0) in OFF Tohoku Area for 1885-1980}
\description{
  Earthquakes data used by Yosihiko Ogata in his 1988 JASA paper.
}
\usage{data(ShallowShocks)}
\format{
  A \code{data.frame} with the following variables:

  \tabular{ll}{
    year: \tab year of occurrence.\cr
    month: \tab month of occurrence.\cr
    day: \tab day of occurrence.\cr
    hour: \tab hour of occurrence.\cr
    minute: \tab minute of occurrence.\cr
    magnitude: \tab magnitude on Richter's scale.\cr
    type: \tab type of earthquake: \code{main} (shock), \code{foreshock},
    \code{aftershock}; according to Utsu.\cr
    Date: \tab date in days starting from January 1st 1885.\cr
    energy.sqrt: \tab square root of the energy expressed in erg.\cr
  }
}
\details{
  Quakes 213 and 214 were given exactly the same dates in Ogata
  (1988). Quake 214 has here been delayed by 1 minute.
}
\source{
  Ogata (1988) Table 1, pp 14-15.
}
\references{
  Ogata, Yosihiko (1988) Statistical Models for Earthquake Occurrences and Residual
  Analysis for Point Processes. \emph{Journal of the American
    Statistical Association} \bold{83}: 9-27.
}
\examples{
data(ShallowShocks)
## Reproduce Fig. 2 of Ogata 1988
layout(matrix(1:3, nrow = 3))
plot(ShallowShocks$Date,
     cumsum(ShallowShocks$energy.sqrt) / 10^13,
     type ="l",
     xlab = "",
     ylab = "",
     main = "Cumulative square root of energy")
plot(ShallowShocks$Date,
     cumsum(1+numeric(dim(ShallowShocks)[1])),
     type ="l",
     xlab = "",
     ylab = "",
     main = "Cumulative number of shocks")
plot(ShallowShocks$Date,
     ShallowShocks$magnitude,
     type = "h",
     ylim = c(5,9),
     xlab = "Time (days)",
     ylab = "",
     main = "Magnitude vs Occurrence time")
}
\keyword{datasets}
