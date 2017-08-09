\name{dllogis}
\alias{dllogis}
\alias{pllogis}
\alias{qllogis}
\alias{rllogis}
\title{The Log Logistic Distribution}
\description{
  Density, distribution function, quantile function,
  and random generation for the log logistic.
}
\usage{
dllogis(x, location = 0, scale = 1, log = FALSE)
pllogis(q, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
qllogis(p, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
rllogis(n, location = 0, scale = 1)
}

\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.} 
  \item{location, scale}{location and scale parameters (non-negative numeric).}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \code{P[X <= x]}, otherwise, \code{P[X > x]}.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities p are given as log(p).}
}
\details{
  If \code{location} or \code{scale} are omitted, they assume the
  default values of 0 and 1 respectively.

  The log-Logistic distribution with \code{location = m} and \code{scale
    = s} has distribution function

  \deqn{\mathrm{F}(x) = \frac{1}{1+ \exp(-\frac{\log (x) - m}{s})}}{F(x) = 1 / (1 + exp(-(log(x)-m)/s))}

  and density

  \deqn{\mathrm{f}(x)=\frac{1}{s \, x} \frac{\exp (-\frac{\log (x) - m}{s})}{(1+ \exp(-\frac{\log (x) - m}{s}))^2}}{f(x) = 1/(s*x) exp(-(log(x)-m)/s) (1 + exp(-(log(x)-m)/s))^-2.}
  
}
\value{
  \code{dllogis} gives the density, \code{pllogis} gives the
  distribution function, \code{qllogis} gives the quantile function
  and \code{rllogis} generates random deviates.
}
\references{
  Lindsey, J.K. (2004) \emph{Introduction to Applied Statistics: A
    Modelling Approach}. OUP.

  Lindsey, J.K. (2004) \emph{The Statistical Analysis of Stochastic
    Processes in Time}. CUP.
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com} }
\seealso{
  \code{\link{llogisMLE}},
  \code{\link{Lognormal}},
  \code{\link{hllogis}}}
\examples{
\dontrun{
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
