\name{coef.durationFit}
\alias{coef.durationFit}
\alias{logLik.durationFit}
\alias{is.durationFit}
\title{Utility Functions for durationFit Objects}
\description{
  \code{coef.durationFit} and \code{logLik.durationFit} extract
  components of a \code{durationFit} object, while \code{is.durationFit}
  tests if its argument is such an object.
}
\usage{
\method{coef}{durationFit}(object,...)
\method{logLik}{durationFit}(object,...)
is.durationFit(obj)
}
\arguments{
  \item{object}{a \code{durationFit} object.}
  \item{obj}{an object to be tested
    against a \code{durationFit} object.}
  \item{\dots}{see \code{\link{coef}} and \code{\link{logLik}}.}
}
\details{
  Everything is trivial here.
}
\value{
  \code{coef.durationFit} returns the coefficients or the estimates or
  the fitted parameters of the object: a 2 elements named vector.

  \code{logLik.durationFit} returns the loglikelihood value.

  \code{is.durationFit} returns \code{TRUE} if its argument is a
  \code{durationFit} object and \code{FALSE} otherwise.
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com} }
\seealso{
  \code{\link{compModels}},
  \code{\link{invgaussMLE}},
  \code{\link{lnormMLE}},
  \code{\link{llogisMLE}},
  \code{\link{rexpMLE}},
  \code{\link{gammaMLE}},
  \code{\link{weibullMLE}}
}
\examples{
\dontrun{
## load CAL1S data
data(CAL1S)
## convert the data into spikeTrain objects
CAL1S <- lapply(CAL1S,as.spikeTrain)
## look at the train of the 1sd neuron
CAL1S[["neuron 1"]]
## fit a invgauss model to the 1st neuron spike train
n1SDFig <- invgaussMLE(CAL1S[["neuron 1"]])
is.durationFit(n1SDFig)
coef(n1SDFig)
logLik(n1SDFig)
}
}
\keyword{distribution}
\keyword{ts}
