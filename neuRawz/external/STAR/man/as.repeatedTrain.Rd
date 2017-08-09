\name{as.repeatedTrain}
\alias{as.repeatedTrain}
\alias{is.repeatedTrain}
\title{Coerce and Test repeatedTrain Objects}
\description{
  \code{as.repeatedTrain} attempts to coerce a list with numeric vector
  elements to a
  \code{repeatedTrain} object while \code{is.repeatedTrain} tests if its
  argument is such an object.
}
\usage{
as.repeatedTrain(x)
is.repeatedTrain(x)
}
\arguments{
  \item{x}{An object to be coerced to or to test against a \code{repeatedTrain} object.}
}
\details{
  A \code{repeatedTrain} object is list of \code{spikeTrain} objects. It
  is used to store the responses of a given neuron to repeated stimulations.
}
\value{
  \code{as.repeatedTrain} returns a \code{repeatedTrain} object or an error.

  \code{is.repeatedTrain} returns \code{TRUE} if its argument is a
  \code{repeatedTrain} object and \code{FALSE} otherwise.
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com} }
\seealso{
  \code{\link{plot.repeatedTrain}},
  \code{\link{print.repeatedTrain}},
  \code{\link{summary.repeatedTrain}},
  \code{\link{psth}},
  \code{\link{raster}},
  \code{\link{as.spikeTrain}},
  \code{\link{is.spikeTrain}}
}
\examples{
## load CAL1V data
data(CAL1V)
## convert them to repeatedTrain objects
CAL1V <- lapply(CAL1V, as.repeatedTrain)
## did the conversion work?
sapply(CAL1V, is.repeatedTrain)
## look at the raster of the 1st neuron
CAL1V[["neuron 1"]]
}
\keyword{ts}
\keyword{survival}
