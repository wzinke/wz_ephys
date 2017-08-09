\name{hgamma}
\alias{hgamma}
\alias{hweibull}
\alias{hlnorm}
\alias{hinvgauss}
\alias{hllogis}
\alias{hrexp}
\title{Hazard Functions for Some Common Duration Distributions}
\description{
  Hazard functions for the gamma, weibull, lognormal, inverse Gaussian,
  log logistic and refractory exponential distributions
}
\usage{
hgamma(x, shape, rate = 1, scale = 1/rate, log = FALSE)
hweibull(x, shape, scale = 1, log = FALSE)
hlnorm(x, meanlog = 0, sdlog = 1, log = FALSE)
hinvgauss(x, mu = 1, sigma2 = 1, boundary = NULL,
          log = FALSE)
hllogis(x, location = 0, scale = 1, log = FALSE)
hrexp(x, rate = 10, rp = 0.005, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of quantiles.}
  \item{shape, scale, rate, sdlog}{strictly positive parameters. See
    corresponding distributions for detail.}
  \item{mu, sigma2, boundary}{parameters associated with the inverse
    Gaussian distribution.}
  \item{meanlog}{parameter associated with the log normal distribution.}
  \item{location, rp}{parameters of the log logistic and refratory exponential.}
  \item{log}{should the log hazard be returned? \code{FALSE} by default.}
}
\details{
  These functions are simply obtained by deviding the density
  by the survival fucntion.
}
\value{
  A vector of hazard rates.
}
\references{
  Lindsey, J.K. (2004) \emph{Introduction to Applied Statistics: A
    Modelling Approach}. OUP.

  Lindsey, J.K. (2004) \emph{The Statistical Analysis of Stochastic
    Processes in Time}. CUP.
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com} }
\seealso{
  \code{\link{dinvgauss}},
  \code{\link{dllogis}},
  \code{\link{drexp}}
}
\examples{
\dontrun{
## use a few plots to compare densities and hazard functions

## lognormal
tSeq <- seq(0.001,0.6,0.001)
meanlog.true <- -2.4
sdlog.true <- 0.4
Yd <- dlnorm(tSeq,meanlog.true,sdlog.true)
Yh <- hlnorm(tSeq,meanlog.true,sdlog.true)
max.Yd <- max(Yd)
max.Yh <- max(Yh)
Yd <- Yd / max.Yd
Yh <- Yh / max.Yh
oldpar <- par(mar=c(5,4,4,4))
plot(tSeq, Yd, type="n", axes=FALSE, ann=FALSE,
     xlim=c(0,0.6), ylim=c(0,1))
axis(2,at=seq(0,1,0.2),labels=round(seq(0,1,0.2)*max.Yd,digits=2))
mtext("Density (1/s)", side=2, line=3)
axis(1,at=pretty(c(0,0.6)))
mtext("Time (s)", side=1, line=3)
axis(4, at=seq(0,1,0.2), labels=round(seq(0,1,0.2)*max.Yh,digits=2))
mtext("Hazard (1/s)", side=4, line=3, col=2)
mtext("Lognormal Density and Hazard Functions", side=3, line=2,cex=1.5)
lines(tSeq,Yd)
lines(tSeq,Yh,col=2)
par(oldpar)

## inverse Gaussian
tSeq <- seq(0.001,0.6,0.001)
mu.true <- 0.075
sigma2.true <- 3
Yd <- dinvgauss(tSeq,mu.true,sigma2.true)
Yh <- hinvgauss(tSeq,mu.true,sigma2.true)
max.Yd <- max(Yd)
max.Yh <- max(Yh)
Yd <- Yd / max.Yd
Yh <- Yh / max.Yh
oldpar <- par(mar=c(5,4,4,4))
plot(tSeq, Yd, type="n", axes=FALSE, ann=FALSE,
     xlim=c(0,0.6), ylim=c(0,1))
axis(2,at=seq(0,1,0.2),labels=round(seq(0,1,0.2)*max.Yd,digits=2))
mtext("Density (1/s)", side=2, line=3)
axis(1,at=pretty(c(0,0.6)))
mtext("Time (s)", side=1, line=3)
axis(4, at=seq(0,1,0.2), labels=round(seq(0,1,0.2)*max.Yh,digits=2))
mtext("Hazard (1/s)", side=4, line=3, col=2)
mtext("Inverse Gaussian Density and Hazard Functions", side=3, line=2,cex=1.5)
lines(tSeq,Yd)
lines(tSeq,Yh,col=2)
par(oldpar)

## gamma
tSeq <- seq(0.001,0.6,0.001)
shape.true <- 6
scale.true <- 0.012
Yd <- dgamma(tSeq, shape=shape.true, scale=scale.true)
Yh <- hgamma(tSeq, shape=shape.true, scale=scale.true)
max.Yd <- max(Yd)
max.Yh <- max(Yh)
Yd <- Yd / max.Yd
Yh <- Yh / max.Yh
oldpar <- par(mar=c(5,4,4,4))
plot(tSeq, Yd, type="n", axes=FALSE, ann=FALSE,
     xlim=c(0,0.6), ylim=c(0,1))
axis(2,at=seq(0,1,0.2),labels=round(seq(0,1,0.2)*max.Yd,digits=2))
mtext("Density (1/s)", side=2, line=3)
axis(1,at=pretty(c(0,0.6)))
mtext("Time (s)", side=1, line=3)
axis(4, at=seq(0,1,0.2), labels=round(seq(0,1,0.2)*max.Yh,digits=2))
mtext("Hazard (1/s)", side=4, line=3, col=2)
mtext("Gamma Density and Hazard Functions", side=3, line=2,cex=1.5)
lines(tSeq,Yd)
lines(tSeq,Yh,col=2)
par(oldpar)

## Weibull
tSeq <- seq(0.001,0.6,0.001)
shape.true <- 2.5
scale.true <- 0.085
Yd <- dweibull(tSeq, shape=shape.true, scale=scale.true)
Yh <- hweibull(tSeq, shape=shape.true, scale=scale.true)
max.Yd <- max(Yd)
max.Yh <- max(Yh)
Yd <- Yd / max.Yd
Yh <- Yh / max.Yh
oldpar <- par(mar=c(5,4,4,4))
plot(tSeq, Yd, type="n", axes=FALSE, ann=FALSE,
     xlim=c(0,0.6), ylim=c(0,1))
axis(2,at=seq(0,1,0.2),labels=round(seq(0,1,0.2)*max.Yd,digits=2))
mtext("Density (1/s)", side=2, line=3)
axis(1,at=pretty(c(0,0.6)))
mtext("Time (s)", side=1, line=3)
axis(4, at=seq(0,1,0.2), labels=round(seq(0,1,0.2)*max.Yh,digits=2))
mtext("Hazard (1/s)", side=4, line=3, col=2)
mtext("Weibull Density and Hazard Functions", side=3, line=2,cex=1.5)
lines(tSeq,Yd)
lines(tSeq,Yh,col=2)
par(oldpar)

## refractory exponential
tSeq <- seq(0.001,0.6,0.001)
rate.true <- 20
rp.true <- 0.01
Yd <- drexp(tSeq, rate.true, rp.true)
Yh <- hrexp(tSeq, rate.true, rp.true)
max.Yd <- max(Yd)
max.Yh <- max(Yh)
Yd <- Yd / max.Yd
Yh <- Yh / max.Yh
oldpar <- par(mar=c(5,4,4,4))
plot(tSeq, Yd, type="n", axes=FALSE, ann=FALSE,
     xlim=c(0,0.6), ylim=c(0,1))
axis(2,at=seq(0,1,0.2),labels=round(seq(0,1,0.2)*max.Yd,digits=2))
mtext("Density (1/s)", side=2, line=3)
axis(1,at=pretty(c(0,0.6)))
mtext("Time (s)", side=1, line=3)
axis(4, at=seq(0,1,0.2), labels=round(seq(0,1,0.2)*max.Yh,digits=2))
mtext("Hazard (1/s)", side=4, line=3, col=2)
mtext("Refractory Exponential Density and Hazard Functions", side=3, line=2,cex=1.5)
lines(tSeq,Yd)
lines(tSeq,Yh,col=2)
par(oldpar)

## log logistic
tSeq <- seq(0.001,0.6,0.001)
location.true <- -2.7
scale.true <- 0.025
Yd <- dllogis(tSeq, location.true, scale.true)
Yh <- hllogis(tSeq, location.true, scale.true)
max.Yd <- max(Yd)
max.Yh <- max(Yh)
Yd <- Yd / max.Yd
Yh <- Yh / max.Yh
oldpar <- par(mar=c(5,4,4,4))
plot(tSeq, Yd, type="n", axes=FALSE, ann=FALSE,
     xlim=c(0,0.6), ylim=c(0,1))
axis(2,at=seq(0,1,0.2),labels=round(seq(0,1,0.2)*max.Yd,digits=2))
mtext("Density (1/s)", side=2, line=3)
axis(1,at=pretty(c(0,0.6)))
mtext("Time (s)", side=1, line=3)
axis(4, at=seq(0,1,0.2), labels=round(seq(0,1,0.2)*max.Yh,digits=2))
mtext("Hazard (1/s)", side=4, line=3, col=2)
mtext("Log Logistic Density and Hazard Functions", side=3, line=2,cex=1.5)
lines(tSeq,Yd)
lines(tSeq,Yh,col=2)
par(oldpar)
}
}
\keyword{distribution}
\keyword{ts}
