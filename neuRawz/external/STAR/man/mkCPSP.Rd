\name{mkCPSP}
\alias{mkCPSP}
\alias{print.CountingProcessSamplePath}
\alias{plot.CountingProcessSamplePath}
\alias{lines.CountingProcessSamplePath}
\alias{as.CPSP}
\title{Counting Process Sample Paths}
\description{
  Functions to create and explore \code{CountingProcessSamplePath}
  objects. These objects are complementary to the \code{spikeTrain}
  objects, the latter being in fact point processes representations.
}
\usage{
mkCPSP(st, from = floor(min(st)), to = ceiling(max(st)))
as.CPSP(x)
\method{print}{CountingProcessSamplePath}(x, digits = 5, ...)
\method{plot}{CountingProcessSamplePath}(x, y, col, lwd, xlim, ylim,
     xlab, ylab, xaxs, yaxs, main, ...)
\method{lines}{CountingProcessSamplePath}(x, ...)
}
\arguments{
  \item{st}{A \code{numeric} vector with \emph{strictly} increasing elements.}
  \item{from}{A \code{numeric}, the time at which the counting process
    obeservation started.}
  \item{to}{A \code{numeric}, the time at which the counting process
    obeservation ended.}
  \item{x}{A \code{numeric} or a \code{spikeTrain} object for
    \code{as.CPSP}, a \code{CountingProcessSamplePath} object for
    \code{print}, \code{plot} and \code{lines}.}
  \item{digits}{An \code{integer}, the number of digits to be used while
  printing summaries. See \code{\link{round}}.}
  \item{y}{Not used but required by the \code{plot} method definition.}
  \item{col,lwd,xlim,ylim,xlab,ylab,main,xaxs,yaxs}{See \code{\link{plot}}.}
  \item{\dots}{Not used in \code{print} (but included for compatibility
    with the method definition) otherwise used like in
    \code{\link{plot}} and \code{\link{lines}}.}
}
\details{
  \code{CountingProcessSamplePath} objects are complementary to
  \code{spikeTrain} objects. They are also used to represente slightly
  more general properties of these objects and are directed towards
  model testing.

  More formaly, if we observe n events at times
  \eqn{\{t_1,\ldots,t_n\}}{{t1,...,tn}} such that, \eqn{from < t_1 <
  \ldots < t_n \le to}{from < t1 < ... < tn <= to}, the \emph{counting
  process sample path} is the right continuous function defined by:
  \deqn{N(t) = \sharp \{ t_j \; \mathrm{with} \; from < t_j \leq t \}}{N(t)
    = # {tj : from < tj <= t}}
  where \eqn{\sharp}{#} stands for the number of elements of a set.

}
\value{
  \code{mkCPSP} returns an object of class
  \code{CountingProcessSamplePath}. This object is a \code{list} with
  the following components:
  \item{cpspFct }{a right continuous \code{function} of \code{t} returning the
    number of events whose occurrence time is strictly larger than
    \code{from} and smaller of equal than \code{t}. \code{t} can be a
    vector. If \code{missing} the cumulative number of events at the
    events occurrence time is returned.}
  \item{ppspFct }{a \code{function} that does not take any argument and
    that returns the sequence of events times, that is, the "point
    process sample path".}
  \item{spikeTrainFct }{a \code{function} that does not take any
    argument and that returns the \code{spikeTrain} object associated
    with the \code{CountingProcessSamplePath} object.}
  \item{from }{argument \code{from} of \code{mkCPSP}.}
  \item{to }{argument \code{to} of \code{mkCPSP}}
  \item{call }{the matched call.}
  
  Functions \code{plot} and \code{lines} are used for their side
  effects, function \code{print} returns a short description of the
  object corresponding to the \code{summary} returned by function
  \code{\link{summary.spikeTrain}} for \code{spikeTrain}
  objects. Function \code{as.CPSP} returns a \code{CountingProcessSamplePath}.
}
\references{
  D. R. Cox and P. A. W. Lewis (1966) \emph{The Statistical Analysis of
    Series of Events}. John Wiley and Sons.
  
  Brillinger, D. R. (1988) Maximum likelihood analysis of spike trains
  of interacting nerve cells. \emph{Biol. Cybern.} \bold{59}: 189--200.

  Johnson, D.H. (1996) Point process models of single-neuron
  discharges. \emph{J. Computational Neuroscience} \bold{3}: 275--299. 
}
\author{ Christophe Pouzat  \email{christophe.pouzat@gmail.com} }
\note{
  This functions are directed towards model testing, don't be surprised
  if they look redundant with the corresponding functions for
  \code{spikeTrain} objects. An apparent difference of detail with the
  latter is that no scale (like seconds) is assumed by default for
  \code{CountingProcessSamplePath} objects. This is to cope in a natural
  way with the time transformation / rescaling procedures used to test
  conditional intensity models. 
}
\seealso{
  \code{\link{summary.CountingProcessSamplePath}},
  \code{\link{print.CountingProcessSamplePath.summary}},
  \code{\link{plot.CountingProcessSamplePath.summary}},
  \code{\link{summary.spikeTrain}},
  \code{\link{print.spikeTrain}},
  \code{\link{plot.spikeTrain}},
  \code{\link{as.spikeTrain}}
  
}
\examples{
## A simple illustration with Ogata's hearthquakes data set
data(ShallowShocks)
plot(mkCPSP(ShallowShocks$Date),
     xlab="Time (days)",
     main="Shallow Shocks Counting Process of Ogata 1988")
## An illustration with on of STAR's data neuroanl dicharge data set
data(e060824spont)
## Create the object from a spikeTrain
n1spt.cp <- as.CPSP(e060824spont[["neuron 1"]])
## print it
n1spt.cp
## plot it
plot(n1spt.cp)
}
\keyword{ts}
\keyword{survival}
