\name{diff.spikeTrain}
\alias{diff.spikeTrain}
\title{ diff method for spikeTrain objects}
\description{
  \code{diff} method for \code{spikeTrain} objects.
}
\usage{
\method{diff}{spikeTrain}(x, ...)
}
\arguments{
  \item{x}{ a \code{spikeTrain} object.}
  \item{\dots}{ see \code{\link{diff}}}
}
\value{
  a \code{numeric}
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com} }
\seealso{
  \code{\link{diff}},
  \code{\link{as.spikeTrain}},
  \code{\link{is.spikeTrain}}
}
\examples{
data(CAL1S)
## convert data into spikeTrain objects
CAL1S <- lapply(CAL1S,as.spikeTrain)
## look at the individual trains
## first the "raw" data
CAL1S[["neuron 1"]]
## get the isi of neuron 1
n1.isi <- diff(CAL1S[["neuron 1"]])
}
\keyword{ts}
