\name{varianceTime}
\alias{varianceTime}
\alias{plot.varianceTime}
\alias{is.varianceTime}
\title{Variance-Time Analysis for Spike Trains}
\description{
  Performs Variance-Time Analysis for a Spike Train (or any univariate
  time series) assuming a Poisson Process with the same Rate as the
  Spike Train.
}
\usage{
varianceTime(spikeTrain, CI = c(0.95, 0.99), windowSizes)
is.varianceTime(obj)
\method{plot}{varianceTime}(x, style = c("default", "Ogata"),
     unit = "s", xlab, ylab, main, sub, xlim, ylim, ...)

}

\arguments{
  \item{spikeTrain}{a \code{spikeTrain} object or a vector which can be
    coerced to such an object.}
  \item{obj}{a object to test against a \code{varianceTime} object.}
  \item{x}{a \code{varianceTime} object.}
  \item{CI}{a numeric vector with at most two elements. The coverage
    probability of the confidence intervals.}
  \item{windowSizes}{a numeric increasing vector of positive
    numbers. The window sizes used to split the spike train.}
  \item{style}{a character. The style of the plot, \code{"default"} or \code{"Ogata"}.}
  \item{unit}{a character. The unit in which the spike times are expressed.}
  \item{xlab}{a character. The x label.}
  \item{ylab}{a character. The y label.}
  \item{main}{a character. The title.}
  \item{sub}{a character. The subtitle.}
  \item{xlim}{a numeric. See \code{\link{plot}}.}
  \item{ylim}{a numeric. See \code{\link{plot}}.}
  \item{\dots}{see \code{\link{plot}}.}
}
\details{
  See Fig. 5 of Ogata (1988) for details. The confidence intervals are
  obtained with a Normal approximation of the Poisson distribution.
}
\value{
  \code{varianceTime} returns a list of class \code{varianceTime} with the following elements:
  
  \item{s2}{numeric vector of empirical variance.}
  \item{sigma2}{numeric vector of expected variance under the Poisson hypothesis.}
  \item{ciUp}{a numeric vector or a 2 rows matrix with the upper limits
    of the confidence interval(s).}
  \item{ciLow}{a numeric vector or a 2 rows matrix with the lower limits
    of the confidence interval(s).}
  \item{windowSizes}{numeric vector of window sizes actually used.}
  \item{CI}{a numeric vector, the coverage probabilities of the confidence intervals.}
  \item{call}{the matched call}

  \code{plot.varianceTime} is used for its side effect: a graph is
  produced.

  \code{is.varianceTime} returns \code{TRUE} if its argument is a
  \code{varianceTime} object and \code{FALSE} otherwise.
}
\references{
  Ogata, Yosihiko (1988) Statistical Models for Earthquake Occurrences and Residual
  Analysis for Point Processes. \emph{Journal of the American
    Statistical Association} \bold{83}: 9-27.
}
\author{
  Christophe Pouzat  \email{christophe.pouzat@gmail.com} and Chong Gu
  \email{chong@stat.purdue.edu} for a correction on the sampling
  variance of the variance of a normal distribution.
}
\seealso{
  \code{\link{acf.spikeTrain}},
  \code{\link{renewalTestPlot}}
}
\examples{
## Replicate (almost) Fig. 5 of Ogata 1988
data(ShallowShocks)
vtShallow <- varianceTime(ShallowShocks$Date,,c(5,10,20,40,60,80,seq(100,500,by = 25))*10)
is.varianceTime(vtShallow)
plot(vtShallow, style="Ogata")
}
\keyword{ts}
\keyword{survival}
