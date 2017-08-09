\name{dinvgauss}
\alias{dinvgauss}
\alias{qinvgauss}
\alias{pinvgauss}
\alias{rinvgauss}
\title{The Inverse Gaussian Distribution}
\description{
  Density, distribution function, quantile function,
  and random generation for the inverse Gaussian.
}
\usage{
dinvgauss(x, mu = 1, sigma2 = 1, boundary = NULL,
          log = FALSE)
pinvgauss(q, mu = 1, sigma2 = 1,
          boundary = NULL, lower.tail = TRUE,
          log.p = FALSE)
qinvgauss(p, mu = 1, sigma2 = 1,
          boundary = NULL)
rinvgauss(n = 1, mu = 1, sigma2 = 1, boundary = NULL)
}

\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is taken to be the number required.}
  \item{mu}{mean value of the distribution in the default
    parameterization, \code{mean value / boundary} otherwise. Can also
    be viewed as the inverse of the drift of the latent Brownian motion.}
  \item{sigma2}{variance of the latent Brownian motion. When this
    parameterization is used (the default) the distance between the "starting" point
    and the boundary ("absorbing barrier") is set to 1.}
  \item{boundary}{distance between the starting point and the "absorbing
    barrier" of the latent Brownian motion. When this parameterization is
    used the Brownian motion variance is set to 1.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \code{P[X <= x]}, otherwise, \code{P[X > x]}.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities p are given as log(p).}
}
\details{
  With the default, \code{"sigma2"}, parameterization (\code{mu = m,
    sigma2 = s^2}) the inverse
  Gaussian distribution has density:
  \deqn{%
    f(x)=\frac{1}{\sqrt{2 \, \pi \, \sigma^2 \, x^3}} \, \exp%
    (-\frac{1}{2}\frac{(x-\mu)^2}{x \, \sigma^2 \, \mu^2})
  }{%
    f(x) = 1/sqrt(2*pi*x^3*s^2) * exp(-0.5*(x-mu)^2/(x*s^2*m^2))
  }
  with \eqn{\sigma^2 > 0}{s^2 > 0}.
  The theoretical mean is: \eqn{\mu}{m} and the theoretical variance is:
  \eqn{\mu^3 \sigma^2}{m^3*s^2}.
  With the default, \code{"boundary"}, parameterization (\code{mu = m,
   boundary = b})the inverse
  Gaussian distribution has density:
  \deqn{%
    f(x)=\frac{b}{\sqrt{2 \, \pi \, x^3}} \, \exp%
    (-\frac{1}{2}\frac{(x-b \, \mu)^2}{x \, \mu^2})
  }{%
    f(x) = (b/sqrt(2*pi*x^3)) * exp(-0.5*(x-b*mu)^2/(x*m^2))
  }
  with \eqn{\sigma^2 > 0}{s^2 > 0}.
  The theoretical mean is: \eqn{\mu \, b}{m * b} and the theoretical variance is:
  \eqn{\mu^3 \sigma^2}{m^3*b}.
  The latent Brownian motion is described in Lindsey (2004) pp 209-213,
  Whitemore and Seshadri (1987), Aalen and Gjessing (2001) and Gerstein
  and Mandelbrot (1964).

  The expression for the distribution function is given in Eq. 4 of
  Whitemore and Seshadri (1987).

  Initial guesses for the inversion of the distribution function used
  in \code{qinvgauss} are obtained with the transformation of Whitemore
  and Yalovsky (1978).
  
  Random variates are obtained with the method of Michael et al (1976)
  which is also described by Devroye (1986, p 148) and Gentle (2003, p 193). 
  
  
}
\value{
  \code{dinvgauss} gives the density, \code{pinvgauss} gives the
  distribution function, \code{qinvgauss} gives the quantile function
  and \code{rinvgauss} generates random deviates.
}
\references{
  Gerstein, George L. and Mandelbrot, Benoit (1964) Random Walk Models
  for the Spike Activity of a Single Neuron. \emph{Biophys J.} \bold{4}:
  41--68. \url{http://www.pubmedcentral.nih.gov/articlerender.fcgi?tool=pubmed&pubmedid=14104072}. 

  Whitemore, G. A. and Yalovsky, M. (1978) A normalizing logarithmic
  transformation for inverse Gaussian random
  variables. \emph{Technometrics} \bold{20}: 207--208.
  
  Whitmore, G. A. and Seshadri, V. (1987) A Heuristic
  Derivation of the Inverse Gaussian Distribution. \emph{The American
    Statistician} \bold{41}: 280--281.

  Aalen, Odd O. and Gjessing, Hakon K. (2001) Understanding the Shape of
  the Hazard Rate: A Process Point of View. \emph{Statistical Science}
  \bold{16}: 1--14.

  Lindsey, J.K. (2004) \emph{Introduction to Applied Statistics: A
    Modelling Approach}. OUP.

  Michael, J. R., Schucany, W. R. and Haas, R. W. (1976) Generating
  random variates using transformations with multiple roots. \emph{The
    American Statistician} \bold{30}: 88--90.

  Devroye, L. (1986) \emph{Non-Uniform Random Variate
    Generation}. Springer-Verlag. \url{http://cg.scs.carleton.ca/~luc/rnbookindex.html}.

  Gentle, J. E. (2003) \emph{Random Number Generation and Monte Carlo
    Methods}. Springer.
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com} }
\seealso{
  \code{\link{invgaussMLE}},
  \code{\link{Lognormal}},
  \code{\link{hinvgauss}}
}
\examples{
\dontrun{
## Start with the inverse Gauss
## Define standard mu and sigma
mu.true <- 0.075 ## a mean ISI of 75 ms
sigma2.true <- 3
## Define a sequence of points on the time axis
X <- seq(0.001,0.3,0.001)
## look at the density
plot(X,dinvgauss(X,mu.true,sigma2.true),type="l",xlab="ISI (s)",ylab="Density")

## Generate a sample of 100 ISI from this distribution
sampleSize <- 100
sampIG <- rinvgauss(sampleSize,mu=mu.true,sigma2=sigma2.true)
## check out the empirical survival function (obtained with the Kaplan-Meyer
## estimator) against the true one
library(survival)
sampIG.KMfit <- survfit(Surv(sampIG,1+numeric(length(sampIG))) ~1)
plot(sampIG.KMfit,log=TRUE)
lines(X,pinvgauss(X,mu.true,sigma2.true,lower.tail=FALSE),col=2)

## Get a ML fit
sampIGmleIG <- invgaussMLE(sampIG)
## compare true and estimated parameters
rbind(est = sampIGmleIG$estimate,se = sampIGmleIG$se,true = c(mu.true,sigma2.true))
## plot contours of the log relative likelihood function
Mu <- seq(sampIGmleIG$estimate[1]-3*sampIGmleIG$se[1],
          sampIGmleIG$estimate[1]+3*sampIGmleIG$se[1],
          sampIGmleIG$se[1]/10)
Sigma2 <- seq(sampIGmleIG$estimate[2]-7*sampIGmleIG$se[2],
              sampIGmleIG$estimate[2]+7*sampIGmleIG$se[2],
              sampIGmleIG$se[2]/10)
sampIGmleIGcontour <- sapply(Mu, function(mu) sapply(Sigma2, function(s2) sampIGmleIG$r(mu,s2))) 
contour(Mu,Sigma2,t(sampIGmleIGcontour),
        levels=c(log(c(0.5,0.1)),-0.5*qchisq(c(0.95,0.99),df=2)),
        labels=c("log(0.5)",
          "log(0.1)",
          "-1/2*P(Chi2=0.95)",
          "-1/2*P(Chi2=0.99)"),
        xlab=expression(mu),ylab=expression(sigma^2))
points(mu.true,sigma2.true,pch=16,col=2)
## We can see that the contours are more parabola like on a log scale
contour(log(Mu),log(Sigma2),t(sampIGmleIGcontour),
        levels=c(log(c(0.5,0.1)),-0.5*qchisq(c(0.95,0.99),df=2)),
        labels=c("log(0.5)",
          "log(0.1)",
          "-1/2*P(Chi2=0.95)",
          "-1/2*P(Chi2=0.99)"),
        xlab=expression(log(mu)),ylab=expression(log(sigma^2)))
points(log(mu.true),log(sigma2.true),pch=16,col=2)
## make a deviance test for the true parameters
pchisq(-2*sampIGmleIG$r(mu.true,sigma2.true),df=2)
## check fit with a QQ plot
qqDuration(sampIGmleIG,log="xy")

## Generate a censored sample using an exponential distribution
sampEXP <- rexp(sampleSize,1/(2*mu.true))
sampIGtime <- pmin(sampIG,sampEXP)
sampIGstatus <- as.numeric(sampIG <= sampEXP)
## fit the censored sample
sampIG2mleIG <- invgaussMLE(sampIGtime,,sampIGstatus)
## look at the results
rbind(est = sampIG2mleIG$estimate,
      se = sampIG2mleIG$se,
      true = c(mu.true,sigma2.true))
pchisq(-2*sampIG2mleIG$r(mu.true,sigma2.true),df=2)
## repeat the survival function estimation
sampIG2.KMfit <- survfit(Surv(sampIGtime,sampIGstatus) ~1)
plot(sampIG2.KMfit,log=TRUE)
lines(X,pinvgauss(X,sampIG2mleIG$estimate[1],sampIG2mleIG$estimate[2],lower.tail=FALSE),col=2)
}
}
\keyword{distribution}
\keyword{ts}
