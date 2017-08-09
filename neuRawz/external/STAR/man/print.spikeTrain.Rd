\name{print.spikeTrain}
\alias{print.spikeTrain}
\alias{summary.spikeTrain}
\title{Print and Summary Methods for spikeTrain Objects}
\description{
  Print and summary \code{\link{methods}} for \code{spikeTrain} objects.
}
\usage{
\method{print}{spikeTrain}(x,...)
\method{summary}{spikeTrain}(object, timeUnit = "s", digits = 3, ...)
}

\arguments{
  \item{x, object}{A \code{spikeTrain} object.}
  \item{timeUnit}{The unit with which the occurrence times were measured.}
  \item{digits}{The number of digits used to print the summary (see
    \code{\link{round}}).}
  \item{\dots}{see \code{\link{print}} and \code{\link{summary}}.}
}
\details{
  \code{print.spikeTrain} does in fact call the \code{plot} method for
  \code{spikeTrain} objects.
}
\value{
  \code{print.spikeTrain} generates a plot as a side effect.

  \code{summary.spikeTrain} returns the number of spikes, the times of
  the first and last spikes, the mean inter-spike interval (ISI) and its sd as
  well as the mean and sd of the log(ISI) together with the shortest and
  longest ISIs.
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com} }
\seealso{
  \code{\link{as.spikeTrain}},
  \code{\link{is.spikeTrain}},
  \code{\link{renewalTestPlot}},
  \code{\link{varianceTime}},
  \code{\link{stepfun}}
}
\examples{
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
\keyword{survival}

