\name{plot.spikeTrain}
\alias{plot.spikeTrain}
\title{Display Counting Process Associated with Single Spike Train}
\description{
  Adds a counting process display to the classical raster plot of single
  spike trains.
}
\usage{
\method{plot}{spikeTrain}(x, xlab = "Time (s)", ylab = "Cumulative Number of Events",
                main = paste("Counting Process of",deparse(substitute(x))),
                xlim = c(floor(x[1]), ceiling(x[length(x)])),
                ylim = c(0, length(x) + 1),
                do.points = ifelse(length(x) < 100, TRUE, FALSE),
                addMeanRate = TRUE, addRug = TRUE, ...)
}
\arguments{
  \item{x}{a \code{spikeTrain} object or a vector which can be
    coerced to such an object.}
  \item{xlab}{a character. The x label.}
  \item{ylab}{a character. The y label.}
  \item{main}{a character. The title.}
  \item{xlim}{a numeric. See \code{\link{plot}}.}
  \item{ylim}{a numeric. See \code{\link{plot}}.}
  \item{do.points}{see \code{\link{plot.stepfun}}.}
  \item{addMeanRate}{should the expected counting process for a Poisson
    process with the same rate be added to the plot?}
  \item{addRug}{should a rug representation be added at teh bottom of
    the plot? See \code{\link{rug}}.}
  \item{\dots}{additional arguments passed to \code{plot}, see
    \code{\link{plot}} and \code{\link{plot.stepfun}}.} 
}
\details{
  The counting process is obtained by a call to \code{\link{stepfun}}.
  When \code{xlab}, \code{ylab}, \code{main}, \code{xlim} or \code{ylim}
  is (are) missing, default values are used.
}
\value{
  Nothing is returned, \code{plot.spikeTrain} is used for its side effect, a
  plot is generated on the current graphic device.
}
\references{
  D. R. Cox and P. A. W. Lewis (1966) \emph{The Statistical Analysis of
    Series of Events}. John Wiley and Sons.
  
  Brillinger, D. R. (1988) Maximum likelihood analysis of spike trains
  of interacting nerve cells. \emph{Biol. Cybern.} \bold{59}: 189--200.

  Johnson, D.H. (1996) Point process models of single-neuron
  discharges. \emph{J. Computational Neuroscience} \bold{3}: 275--299. 
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com} }
\seealso{
  \code{\link{as.spikeTrain}},
  \code{\link{is.spikeTrain}},
  \code{\link{print.spikeTrain}},
  \code{\link{summary.spikeTrain}},
  \code{\link{renewalTestPlot}},
  \code{\link{varianceTime}},
  \code{\link{stepfun}},
  \code{\link{plot.stepfun}},
  \code{\link{rug}}
}
\examples{
\dontrun{
data(ShallowShocks)
plot(as.spikeTrain(ShallowShocks$Date),
     xlab="Time (days)",
     main="Shallow Shocks Counting Process of Ogata 1988")
}
}
\keyword{ts}
\keyword{survival}
