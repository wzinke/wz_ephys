\name{drexp}
\alias{drexp}
\alias{prexp}
\alias{qrexp}
\alias{rrexp}
\title{The Refractory Exponential Distribution}
\description{
  Density, distribution function, quantile function,
  and random generation for the refractory exponential.
}
\usage{
drexp(x, rate = 10, rp = 0.005, log = FALSE)
prexp(q, rate = 10, rp = 0.005, lower.tail = TRUE, log.p = FALSE)
qrexp(p, rate = 10, rp = 0.005, lower.tail = TRUE, log.p = FALSE)
rrexp(n, rate = 10, rp = 0.005)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.} 
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \code{P[X <= x]}, otherwise, \code{P[X > x]}.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities p are given as log(p).}
  \item{rate}{rate parameter (non-negative numeric).}
  \item{rp}{refractory period parameter (non-negative numeric).}
}
\details{
  The refractory exponential distribution with \code{rate}, r, and \code{refractory
  period}, rp, has density:

  f(x) = r exp(- r (x-rp))

  for \code{x >= rp}. 
}
\value{
  \code{drexp} gives the density, \code{prexp} gives the
  distribution function, \code{qrexp} gives the quantile function
  and \code{rrexp} generates random deviates.
}
\references{
  Johnson, D. H. and Swami, A. (1983) The transmission of signals by
  auditory-nerve fiber discharge patterns. \emph{J. Acoust. Soc. Am.}
  \bold{74}: 493--501.
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com} }
\seealso{
  \code{\link{rexpMLE}}
}
\examples{
\dontrun{
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
}
}
\keyword{distribution}
\keyword{ts}
