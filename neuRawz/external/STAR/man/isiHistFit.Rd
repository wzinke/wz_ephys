\name{isiHistFit}
\alias{isiHistFit}
\title{ISI Histogram With Fitted Model and CI}
\description{
  Fits a duration model to isis from a spike train. Confidence intervals
  are also drawn.
}
\usage{
isiHistFit(spikeTrain, model, nbins = 10, CI = 0.95, ...)
}
\arguments{
  \item{spikeTrain}{a \code{spikeTrain} object or a numeric vector that
    can be coerced to such an object.}
  \item{model}{a character vector whose elements are selected among:
    \code{"invgauss"}, \code{"lnorm"}, \code{"gamma"},
    \code{"weibull"}, \code{"llogis"}, \code{"rexp"}.}
  \item{nbins}{the number of bins to use.}
  \item{CI}{the confidence coefficient.}
  \item{\dots}{additional arguments passed to \code{hist}, see
    \code{\link{hist}}.}
}
\details{
  Assuming that the train is reasonably well described by a renewal process,
  a \code{model} distribution is fitted to the inter-spike intervals
  (isis) obtained from \code{spikeTrain}. The fitted distribution is
  then used to set the histogram breaks such that a uniform bin
  count would be expected if the fitted distribution was the true
  one. Confidence segments are also obtained from the binomial
  distribution. The histogram is build and the fitted density together
  with confidence intervals are drawn.
}
\value{
  Nothing returned, \code{isiHistFit} is used for its side effect, a
  plot is generated on the current graphic device.
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com} }
\seealso{
  \code{\link{compModels}},
  \code{\link{hist}}
}
\examples{
\dontrun{
## load spontaneous data of 4 putative projection neurons
## simultaneously recorded from the cockroach (Periplaneta
## americana) antennal lobe
data(CAL1S)
## convert data into spikeTrain objects
CAL1S <- lapply(CAL1S,as.spikeTrain)
## look at the individual trains
## first the "raw" data
CAL1S[["neuron 1"]]
## next some summary information
summary(CAL1S[["neuron 1"]])
## next the renewal tests
renewalTestPlot(CAL1S[["neuron 1"]])
## It does not look too bad so let fit simple models
compModels(CAL1S[["neuron 1"]])
## the best one is the invgauss. Let's look at
## it in detail
isiHistFit(CAL1S[["neuron 1"]],"invgauss",xlim=c(0,0.5))
}
}
\keyword{distribution}
\keyword{ts}
