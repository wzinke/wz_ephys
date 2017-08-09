\name{gamlockedTrain}
\alias{gamlockedTrain}
\alias{print.gamlockedTrain}
\alias{summary.gamlockedTrain}
\alias{plot.gamlockedTrain}
\title{Function to Smooth a lockedTrain Object and Related Methods: The
  Penalized Regression Spline Approach}
\description{
  Smooths a \code{lockedTrain} object using a \code{gam} model with the Poisson
  family after binning the object.
}
\usage{
gamlockedTrain(lockedTrain, bw = 0.001, bs = "cr", k = 100, ...)
\method{print}{gamlockedTrain}(x, ...)
\method{summary}{gamlockedTrain}(object, ...)
\method{plot}{gamlockedTrain}(x, xlab, ylab, main, xlim, ylim, col, lwd, ...)

}
\arguments{
  \item{lockedTrain}{a \code{\link{lockedTrain}} object.}
  \item{bw}{the bin width (in s) used to generate the observations on which
  the gam fit will be performed. See details below.}
  \item{bs}{the type of splines used. See \code{\link[mgcv]{s}}.}
  \item{k}{the dimension of the basis used to represent the smooth
    psth. See \code{\link[mgcv]{s}}.}
  \item{x}{an \code{gamlockedTrain} object.}
  \item{object}{an \code{gamlockedTrain} object.}
  \item{xlim}{a numeric (default value supplied). See
    \code{\link{plot}}.}
  \item{ylim}{a numeric (default value supplied). See \code{\link{plot}}.}
  \item{xlab}{a character (default value supplied). See \code{\link{plot}}.}
  \item{ylab}{a character (default value supplied). See \code{\link{plot}}.}
  \item{main}{a character (default value supplied). See \code{\link{plot}}.}
  \item{lwd}{line width used to plot the estimated density. See \code{\link{plot}}.}
  \item{col}{color used to plot the estimated density. See \code{\link{plot}}.}
  \item{\dots}{additional arguments passed to \code{\link[mgcv]{gam}} in
  \code{gamlockedTrain}. Not used in \code{print.gamlockedTrain} and
  \code{summary.gamlockedTrain}. Passed to \code{\link{plot}} in
  \code{plot.gamlockedTrain}.}
}
\details{
  \code{gamlockedTrain} essentially generates a smooth version of the
  histogram obtained by \code{\link{hist.lockedTrain}}. The Idea is to
  build the histogram first with a "too" small bin width before fitting
  a regression spline to it with a Poisson distribution of the observed
  counts.
}
\value{
  A list of class \code{gamlockedTrain} is returned by
  \code{gamlockedTrain}. This list has the following components:
  \item{gamFit }{the \code{\link[mgcv]{gamObject}} generated.}
  \item{Time }{the vector of bin centers.}
  \item{nRef }{the number of spikes in the reference train. See
    \code{\link{hist.lockedTrain}}.}
  \item{testFreq }{the mean frequency of the test neuron. See
    \code{\link{hist.lockedTrain}}.}
  \item{bwV }{the vector of bin widths used.}
  \item{CCH }{a logical which is \code{TRUE} if a cross-intensity was
    estimated and \code{FALSE} in the case of an auto-intensity.}
  \item{call }{the matched call.}

  \code{print.gamlockedTrain} returns the result of \code{\link[mgcv]{print.gam}}
  applied to the component \code{gamFit} of its argument.

  \code{summary.gamlockedTrain} returns the result of \code{\link[mgcv]{summary.gam}}
  applied to the component \code{gamFit} of its argument.

}
\references{
  Wood S.N. (2006) \emph{Generalized Additive Models: An Introduction
    with R}. Chapman and Hall/CRC Press.
}
\author{Christophe Pouzat \email{christophe.pouzat@gmail.com} }
\seealso{
  \code{\link{lockedTrain}},
  \code{\link{plot.lockedTrain}},
  \code{\link[mgcv]{gam}}
}
\examples{
\dontrun{
## load e070528spont data set
data(e070528spont)
## create a lockedTrain object with neuron 1 as reference
## and neuron 3 as test up to lags of +/- 250 ms
lt1.3 <- lockedTrain(e070528spont[[1]],e070528spont[[3]],laglim=c(-1,1)*0.25)
## look at the cross raster plot
lt1.3
## build a histogram of it using a 10 ms bin width
hist(lt1.3,bw=0.01)
## do it the smooth way
slt1.3 <- gamlockedTrain(lt1.3)
plot(slt1.3)
## do some check on the gam fit
summary(slt1.3)
gam.check(gamObj(slt1.3))
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{smooth}
\keyword{regression}
