\name{acf.spikeTrain}
\alias{acf.spikeTrain}
\title{Auto- Covariance and -Correlation Function Estimation
  for Spike Train ISIs}
\description{
  The function \code{acf.spikeTrain} computes (and by default plots) estimates of the
  autocovariance or autocorrelation function of the inter-spike
  intervals of a spike train. 
}
\usage{
acf.spikeTrain(spikeTrain, lag.max = NULL,
       type = c("correlation", "covariance", "partial"),
       plot = TRUE, na.action = na.fail,
       demean = TRUE, xlab = "Lag (in isi #)",
       ylab = "ISI acf",
       main, ...)
}
\arguments{
  \item{spikeTrain}{a \code{spikeTrain} object or a vector which can be
    coerced to such an object.}
  \item{lag.max}{maximum lag at which to calculate the acf.
    Default is \eqn{10\log_{10}(N)}{10*log10(N)} where \eqn{N} is the
    number of ISIs.  Will
    be automatically limited to one less than the number of ISIs
    in the spike train.}
  \item{type}{character string giving the type of acf to be computed.
    Allowed values are
    \code{"correlation"} (the default), \code{"covariance"} or
    \code{"partial"}.}
  \item{plot}{logical. If \code{TRUE} (the default) the acf is plotted.}
  \item{na.action}{function to be called to handle missing
    values. \code{na.pass} can be used.}
  \item{demean}{logical.  Should the covariances be about the sample
    means?}
  \item{xlab}{x axis label.}
  \item{ylab}{y axis label.}
  \item{main}{title for the plot.}
  \item{\dots}{further arguments to be passed to \code{plot.acf}.}
}
\details{
  Just a wrapper for \code{\link{acf}} function. The first argument,
  \code{spikeTrain}, is processed first to extract the inter-spike
  intervals. \code{acf.spikeTrain} is mainly used to plot what Perkel et
  al (1967) call the \emph{serial correlation coefficient} (Eq. 8) or
  \emph{serial covariance coefficient} (Eq. 7), p 400.
}
\value{
  An object of class \code{"acf"}, which is a list with the following
  elements:

  \item{lag}{A three dimensional array containing the lags at which
    the acf is estimated.}
  \item{acf}{An array with the same dimensions as \code{lag} containing
    the estimated acf.}
  \item{type}{The type of correlation (same as the \code{type}
    argument).}
  \item{n.used}{The number of observations in the time series.}
  \item{series}{The name of the series \code{x}.}
  \item{snames}{The series names for a multivariate time series.}
  
  The lag \code{k} value returned by \code{ccf(x,y)} estimates the
  correlation between \code{x[t+k]} and \code{y[t]}.

  The result is returned invisibly if \code{plot} is \code{TRUE}.
}
\references{
  Perkel D. H., Gerstein, G. L. and Moore G. P. (1967) Neural
  Spike Trains and Stochastic Point Processes. I. The Single Spike
  Train. \emph{Biophys. J.}, \bold{7}:
  391-418. \url{http://www.pubmedcentral.nih.gov/articlerender.fcgi?tool=pubmed&pubmedid=4292791}
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com}}
\seealso{\code{\link{acf}},
  \code{\link{varianceTime}},
  \code{\link{renewalTestPlot}}
}
\examples{
## Simulate a log normal train
train1 <- c(cumsum(rlnorm(301,log(0.01),0.25)))
train1 <- as.spikeTrain(train1)
## Get its isi acf
acf.spikeTrain(train1,lag.max=100)
}
\keyword{ts}
