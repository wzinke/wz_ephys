\name{hist.lockedTrain}
\alias{plot.hist.lockedTrain}
\alias{hist.lockedTrain}
\title{Auto- and Cross-Intensity Function Estimate for Spike Trains}
\description{
  \code{hist.lockedTrain} constructs and \code{plot.hist.lockedTrain}
  plots estimates of what Cox and Lewis (1966) call the auto- or
  cross-intensity functions. The auto-intensity function is also called
  the renewal density by Cox and Lewis and by Perkel et al (1967) while
  it is called the intensity of the Palm distribution by Ogata
  (1988). The (estimate of) the cross-intensity function is called
  cross-correlation function by Perkel et al (1967b) and
  cross-correlation histogram by Brillinger et al (1976). 
}
\usage{
\method{hist}{lockedTrain}(x, bw, breaks = NULL, plot = TRUE, ...)
\method{plot}{hist.lockedTrain}(x, style = c("Ogata", "Brillinger"),
                 CI, unit = "s", xlab, ylab, xlim, ylim,
                 type, pch, ...)
}

\arguments{
  \item{x}{a \code{lockedTrain} object returned by the
    \code{\link{lockedTrain}} function.}
  \item{bw}{a non-negative numeric, the bin width.}
  \item{breaks}{a vector giving the breakpoints between cells. If
    \code{NULL} (default) breaks are built using argument
    \code{bw} and component \code{laglim} of \code{obj}.}
  \item{plot}{a logical. If \code{TRUE} a plot is generated as a side
    effect and nothing is returned, if \code{FALSE} a list of class
    \code{hist.lockedTime} is returned.}
  \item{style}{a character. The style of the plot, \code{"Ogata"} or \code{"Brillinger"}.}
  \item{CI}{a numeric vector with at most two elements. The coverage
    probability of the confidence intervals.}
  \item{unit}{a character. The unit in which the spike times are expressed.}
  \item{xlab}{a character. The x label. Default supplied.}
  \item{ylab}{a character. The y label. Default supplied.}
  \item{xlim}{a numeric. See \code{\link{plot}}. Default supplied.}
  \item{ylim}{a numeric. See \code{\link{plot}}. Default supplied.}
  \item{type}{see \code{\link{lines}}. Default supplied.}
  \item{pch}{see \code{\link{plot}}. Default supplied.}
  \item{\dots}{see \code{\link{plot}}.}
}
\details{
  The intensity of the Palm distribution (Ogata, 1988, p 13) is
  estimated by:
  \deqn{\mathrm{m}(s) = \frac{\mathrm{Prob}(\mathrm{event \, in}
      \, (t+s,t+s+\Delta s) \, \mid \, \mathrm{event \, at} \, t)}{\Delta
      s}}{m(s) = Prob(event in (t+s,t+s+bw) | event at t) / bw}
  It is called \emph{renewal density} by Perkel et al (1967) and defined by
  their Eq. 10, p 404.
  Under the null hypothesis of a stationary Poisson process it is a constant whose value
  is the mean discharge rate.

  The cross-intensity function of two spike trains A and B is estimated by
  (Perkel et al, 1967b, p424, Eq. 4 and 5):
  \deqn{\mathrm{m}_{AB}(s) = \frac{\mathrm{Prob}(\mathrm{B \, event \, in}
      \, (t+s,t+s+\Delta s) \, \mid \, \mathrm{A \, event \, at} \, t)}{\Delta
      s}}{m(s;AB) = Prob(B event in (t+s,t+s+bw) | A event at t) / bw}
  
  The \code{style} argument of \code{plot.hist.lockedTrain} generates a
  plot looking like Fig. 6, p 18 of Ogata (1988) if set to
  \code{"Ogata"}. Using \code{style} \code{"Brillinger"} plots the
  square root of the estimate.
}
\value{
  When argument \code{plot} in \code{hist.lockedTrain} is set to \code{FALSE}
  a list of class \code{hist.lockedTrain} with the following components is
  returned:
  
  \item{density}{the density estimate. Equivalent of the component
    \code{density} returned by \code{\link{hist}}.}
  \item{breaks}{a numeric vector with the breaks in between which spikes
  were counted. Similar to the component of the same name returned by \code{\link{hist}}.}
  \item{mids}{a numeric vector with the mid points of
    \code{breaks}. . Similar to the component of the same name returned
    by \code{\link{hist}}.}
  \item{bw}{the bin width used.}
  \item{nRef }{the total number of reference spikes used.}
  \item{refFreq }{the mean frequency of the reference neuron.}
  \item{testFreq }{the mean frequency of the test neuron.}
  \item{obsTime }{the total observation time used (in s).}
  \item{CCH }{a logical which is \code{TRUE} if a cross-intensity was
    estimated and \code{FALSE} in the case of an auto-intensity.}
  \item{call}{the matched call.}
}
\references{
  Ogata, Yosihiko (1988) Statistical Models for Earthquake Occurrences and Residual
  Analysis for Point Processes. \emph{Journal of the American
    Statistical Association} \bold{83}: 9-27.

  D. R. Cox and P. A. W. Lewis (1966) \emph{The Statistical Analysis of
    Series of Events}. John Wiley and Sons.

  J. A. McFadden (1962) On the Lengths of Intervals in a Stationary
  Point Process. \emph{Journal of the Royal Statistical Society. Series
    B}, \bold{24}: 364-382
  
  Perkel D. H., Gerstein, G. L. and Moore G. P. (1967) Neural
  Spike Trains and Stochastic Point Processes. I. The Single Spike
  Train. \emph{Biophys. J.}, \bold{7}:
  391-418. \url{http://www.pubmedcentral.nih.gov/articlerender.fcgi?tool=pubmed&pubmedid=4292791}

  Perkel D. H., Gerstein, G. L. and Moore G. P. (1967b) Neural
  Spike Trains and Stochastic Point Processes. I. Simultaneous Spike
  Trains. \emph{Biophys. J.}, \bold{7}:
  419-440. \url{http://www.pubmedcentral.nih.gov/articlerender.fcgi?tool=pubmed&pubmedid=4292792}

  David R. Brillinger, Hugh L. Bryant and Jose P. Segundo (1976)
  Identification of synaptic interactions. \emph{Biol Cybern},
  \bold{22}: 213-228.

  David R. Brillinger (1976) Estimation of the Second-Order Intensities
  of a Bivariate Stationary Point Process. \emph{Journal of the Royal
    Statistical Society. Series B (Methodological)}, \bold{38}, 60-66.
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com} }
\note{
  The confidence intervals are obtained from a Poisson distribution with
  parameter:
  \code{refFreq} * \code{testFreq} * \code{bw} * \code{obsTime}. Once the
  quantiles of the Poisson distribution have been obtained they are divided
  by: \code{refFreq} * \code{bw} * \code{obsTime}

  These intervals are valid under the stationary Poisson null hypothesis
  for the auto-intensity estimates. They are valid under the
  \emph{stationary independent} null hypothesis for the
  cross-intensity. \emph{There is NO NEED to assume that the test train
    is Poisson or renewal.} See Perkel et al (1967b) and McFadden (1962)
  for a justification/proof of that. The square root transform of
  Brillinger (1976) and Brillinger et al (1976) is (in my opinion) a
  perfect example of shooting at a sparrow with a bazooka. An oversized
  method to get at an intuitively obvious result. There is moreover
  no need to stabilize the variance if your testing against a Poisson
  with a constant rate since then the variance of the null
  hypothesis is stable to start with. These (square root) transforms are
  useful for least square fits with a Poisson noise but NOT in the present context.
}
\seealso{
  \code{\link{varianceTime}},
  \code{\link{renewalTestPlot}},
  \code{\link{lockedTrain}}
}
\examples{
## Reproduce Fig. 6 of Ogata 1988
data(ShallowShocks)
shalShocks <- lockedTrain(as.spikeTrain(ShallowShocks$Date),,c(0,500))
shalShocksH <- hist(shalShocks,5,plot=FALSE)
plot(shalShocksH,"Ogata",c(0.95,0.99),xlab="TIME LAG (DAYS)",ylab="NUMBER OF EVENTS PER DAY")
## Reproduce Fig. 7 of Ogata 1988
myBinNb <- 101
myMidPoints <- seq(from = 1, to = 6, length.out=myBinNb)
myMidPoints <- 10^myMidPoints/200
myBreaks <- c(0,myMidPoints[-myBinNb] + diff(myMidPoints) / 2)
shalShocksH2 <- hist(shalShocks,breaks=myBreaks,plot=FALSE)
yy <- abs(shalShocksH2$density - shalShocksH2$refFreq)
plot(shalShocksH2$mids[shalShocksH2$density>0],
     yy[shalShocksH2$density>0],
     pch = 1,
     xlim = c(0.001,10000),
     log = "xy",
     xlab = "TIME LAG (DAYS)",
     ylab = "NUMBER OF EVENTS PER DAY"
     )

}
\keyword{ts}
\keyword{survival}
