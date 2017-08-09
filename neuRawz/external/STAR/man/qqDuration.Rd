\name{qqDuration}
\alias{qqDuration}
\title{Quantile-Quantile Plot For Fitted Duration Distributions}
\description{
  Produces a QQ plot of empirical against theoretical
  quantiles of one of the following duration distributions: inverse
  Gaussian, log normal, log logistic, refractory exponential, gamma, weibull.
}
\usage{
qqDuration(durationFit, CI = c(0.95, 0.99),
           type = "l", xlab, ylab, main, sub,
           ylim, dataLwd = 2, ablineCol = 2, ...)
}
\arguments{
  \item{durationFit}{a \code{durationFit} object, that is, a list returned by one of these functions:
    \code{\link{invgaussMLE}}, \code{\link{lnormMLE}}, 
    \code{\link{llogisMLE}}, \code{\link{rexpMLE}},
    \code{\link{gammaMLE}}, \code{\link{weibullMLE}}.}
  \item{CI}{a numeric vector with at most tow components, the confidence
  intervals to be drawn. If \code{NULL}, intervals are not drawn.}
  \item{type, xlab, ylab, main, sub, ylim}{see \code{\link{plot}},
    default values are provided if arguments are missing.}
  \item{dataLwd}{non negative integer, the width of the line used to
    draw the data.}
  \item{ablineCol}{color of the diagonal.}
  \item{\dots}{additional arguments passed to \code{\link{plot}}.}
}
\details{
  If the data to which the model was fitted have censored events, the
  latter are not used to build the empirical quantiles.
}
\value{
  Nothing is returned, the function is used for its side effect, a plot
  is generated.
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com} }
\seealso{
  \code{\link{compModels}},
  \code{\link{invgaussMLE}},
  \code{\link{lnormMLE}},
  \code{\link{llogisMLE}},
  \code{\link{rexpMLE}},
  \code{\link{gammaMLE}},
  \code{\link{weibullMLE}}
}
\examples{
\dontrun{
## Simulate a sample with 100 events from an inverse Gaussian
set.seed(1102006,"Mersenne-Twister")
mu.true <- 0.075
sigma2.true <- 3
sampleSize <- 100
sampIG <- rinvgauss(sampleSize,mu=mu.true,sigma2=sigma2.true)
## Fit it with an inverse Gaussian Model
sampIGmleIG <- invgaussMLE(sampIG)
## draw the QQ plot on a log scale
qqDuration(sampIGmleIG,log="xy")
## Fit it with a log normal Model
sampIGmleLN <- lnormMLE(sampIG)
## draw the QQ plot on a log scale
qqDuration(sampIGmleLN,log="xy")
## Fit it with a gamma Model
sampIGmleGA <- gammaMLE(sampIG)
## draw the QQ plot on a log scale
qqDuration(sampIGmleGA,log="xy")
## Fit it with a Weibull Model
sampIGmleWB <- weibullMLE(sampIG)
## draw the QQ plot on a log scale
qqDuration(sampIGmleWB,log="xy")
## Fit it with a refractory exponential Model
sampIGmleRE <- rexpMLE(sampIG)
## draw the QQ plot on a log scale
qqDuration(sampIGmleRE,log="xy")
## Fit it with a log logisitc Model
sampIGmleLL <- llogisMLE(sampIG)
## draw the QQ plot on a log scale
qqDuration(sampIGmleLL,log="xy")
}
}
\keyword{distribution}
\keyword{ts}
