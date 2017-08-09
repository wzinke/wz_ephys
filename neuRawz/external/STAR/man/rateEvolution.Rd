\name{rateEvolution}
\alias{rateEvolution}
\alias{plot.rateEvolution}
\title{Evaluates and Plots a Spike Train Firing Rate's Evolution}
\description{
  \code{rateEvolution} evaluates and \code{plot.rateEvolution} plots the
  firing rate evolution of a \code{spikeTrain} object. The evaluation is
  done by convolving the spike train with a kernel like in
  \code{density} estimation.
}
\usage{
rateEvolution(x, bw, kernel = c("gaussian", "epanechnikov",
                                "rectangular", "triangular",
                                "biweight", "cosine", "optcosine"),
              n = 512, from, to, na.rm = FALSE, ...)
\method{plot}{rateEvolution}(x, main = NULL, xlab = NULL, ylab = "Rate (Hz)",
                   type = "l", zero.line = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a \code{spikeTrain} object or an object which can be coerced
    to it for \code{rateEvolution} or a
    \code{rateEvolution} object for \code{plot.rateEvolution}.}
  \item{bw}{the kernel bin width in seconds. If \code{missing} it is set
  to 10 times the median inter-spike interval of \code{x}.}
  \item{kernel}{see \code{\link{density}}.}
  \item{n}{see \code{\link{density}}.}
  \item{from}{see \code{\link{density}}.}
  \item{to}{see \code{\link{density}}.}
  \item{na.rm}{see \code{\link{density}}.}
  \item{main}{see \code{plot.density}.}
  \item{xlab}{see \code{plot.density}.}
  \item{ylab}{see \code{plot.density}.}
  \item{type}{see \code{plot.density}.}
  \item{zero.line}{see \code{plot.density}.}
  \item{\dots}{see \code{\link{density}} and \code{plot.density}.}
}
\details{
  \code{rateEvolution} is mainly a wrapper for \code{\link{density}}
  which also adjusts the result of the latter such that the y component
  of the returned list is an instantaneous firing rate. If the length of
  \code{x} is smaller or equal to 1 and if \code{from} or
  \code{to} is (are) \code{missing} the returned object has then each of its
  components set to \code{NA} except \code{data.name} (see below). If
  the length of \code{x} is smaller or equal to 1 and if both
  \code{from} and \code{to} are specified a \code{missing} \code{bw} is
  then set to 3 times the spacing between the points of the regular grid
  on which the density is evaluated.

  \code{plot.rateEvolution} is also a wrapper for
  \code{\link{plot.density}} which only adjust the default value of some arguments.
}
\value{
  \code{rateEvolution} returns a LIST of class \code{rateEvolution}
  which inherits from class \code{density}.
  \item{x}{the \code{n} coordinates of the points where the density is
    estimated. See \code{\link{density}}.}
  \item{y}{the estimated rate (in 1/s).  These will be non-negative,
    but can be zero.}
  \item{bw}{the bandwidth used.}
  \item{n}{the sample size after elimination of missing values.}
  \item{call}{the call which produced the result.}
  \item{data.name}{the deparsed name of the \code{x} argument.}
  \item{has.na}{logical, for compatibility (always \code{FALSE}).}

  \code{plot.rateEvolution} is called for its side effect: a plot is generated.
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com} }
\seealso{
  \code{\link{as.spikeTrain}},
  \code{\link{density}},
  \code{\link{plot.density}},
  \code{\link{mkREdf}}
}
\examples{
## load Purkinje cell data recorded in cell-attached mode
data(sPK)
## coerce sPK to a spikeTrain object
sPK <- lapply(sPK, as.spikeTrain)
## get the rate evolution in ctl condition
sPKreCTL <- rateEvolution(sPK[["ctl"]])
## plot the result
plot(sPKreCTL)
## check the bin width which was actually used
sPKreCTL$bw
## look at the effect of a 10 times larger bw
plot(rateEvolution(sPK[["ctl"]],bw=10*sPKreCTL$bw))
## look at the effect of a 10 times smaller one
plot(rateEvolution(sPK[["ctl"]],bw=sPKreCTL$bw/10))
## get the rate evolution in bicuculline conditions
sPKreBICU <- rateEvolution(sPK[["bicu"]])
## plot results
plot(sPKreBICU,col=2)
## add the ctl rate evolution
lines(sPKreCTL)
}
\keyword{ts}
