\name{as.spikeTrain}
\alias{as.spikeTrain}
\alias{is.spikeTrain}
\alias{[.spikeTrain}
\title{Coerce, Test and Extract from spikeTrain Objects }
\description{
  \code{as.spikeTrain} attempts to coerce a numeric vector to a
  \code{spikeTrain} object while \code{is.spikeTrain} tests if its
  argument is such an object. \code{[.spikeTrain}, extracts a subset of
  a \code{spikeTrain} object.
}
\usage{
as.spikeTrain(x)
is.spikeTrain(x)
\method{[}{spikeTrain}(x,i)
}
\arguments{
  \item{x}{An object to be coerced to or to test against a
  \code{spikeTrain} object or a \code{spikeTrain} object for \code{[}.}
  \item{i}{indices specifying elements to extract. \emph{No gaps are allowed}.}
}
\details{
  A \code{spikeTrain} object is a \code{numeric} vector whose elements
  are strictly increasing (that is, something which can be interpreted
  as a sequence of times of successive events with no two events
  occurring at the same time). The extractor method, \code{[} requires
that the extracted elements are without gaps, an error is returned otherwise.
}
\value{
  \code{as.spikeTrain} returns a \code{spikeTrain} object or an error.

  \code{is.spikeTrain} returns \code{TRUE} if its argument is a
  \code{spikeTrain} object and \code{FALSE} otherwise.

  \code{[} returns a \code{spikeTrain} object or an error.
}
\references{
  Perkel D. H., Gerstein, G. L. and Moore G. P. (1967) Neural
  Spike Trains and Stochastic Point Processes. I. The Single Spike
  Train. \emph{Biophys. J.}, \bold{7}:
  391-418. \url{http://www.pubmedcentral.nih.gov/articlerender.fcgi?tool=pubmed&pubmedid=4292791}
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com} }
\seealso{
  \code{\link{plot.spikeTrain}},
  \code{\link{print.spikeTrain}},
  \code{\link{summary.spikeTrain}}
}
\examples{
## load CAL1S data
data(CAL1S)
## convert the data into spikeTrain objects
CAL1S <- lapply(CAL1S,as.spikeTrain)
## Are the list eleemnts now spikeTrain objects?
sapply(CAL1S, is.spikeTrain)
## look at the train of the 1st neuron
CAL1S[["neuron 1"]]
## look at the window 10-40 using the extractor function
CAL1S[["neuron 1"]][10 < CAL1S[["neuron 1"]] & CAL1S[["neuron 1"]] < 40] 

}
\keyword{ts}
\keyword{survival}

