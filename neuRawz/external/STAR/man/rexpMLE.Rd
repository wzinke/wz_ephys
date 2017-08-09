\name{rexpMLE}
\alias{rexpMLE}
\title{Maximum Likelihood Parameter Estimation of a Refractory Exponential Model with Possibly
  Censored Data}
\description{
  Estimate refractory exponential model parameters by the maximum likelihood
  method using possibly censored data.
}
\usage{
rexpMLE(yi, ni = numeric(length(yi)) + 1,
        si = numeric(length(yi)) + 1)
}
\arguments{
  \item{yi}{vector of (possibly binned) observations or a
    \code{spikeTrain} object.}
  \item{ni}{vector of counts for each value of \code{yi}; default: \code{numeric(length(yi))+1}.}
  \item{si}{vector of counts of \emph{uncensored} observations for each
    value of \code{yi}; default: \code{numeric(length(yi))+1}.}
}
\details{
  The MLE are available in closed form even in the censored case for
  this model. The likelihood function cannot be differentiated with
  respect to the \code{rp} (refractory period) parameter at the
  maximum. COnfidence intervals for this parameter are therefore not
  available.

}
\value{
  A list of class \code{durationFit} with the following components:
  \item{estimate}{the estimated parameters, a named vector.}
  \item{se}{the standard errors, a named vector.}
  \item{logLik}{the log likelihood at maximum.}
  \item{r}{a function returning the log of the relative likelihood function.}
  \item{mll}{a function returning the opposite of the log likelihood
    function using the log of the parameters.}
  \item{call}{the matched call.}
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com} }
\seealso{
  \code{\link{drexp}},
  \code{\link{invgaussMLE}},
  \code{\link{lnormMLE}},
  \code{\link{gammaMLE}},
  \code{\link{weibullMLE}}
}
\examples{
\dontrun{
## Simulate sample of size 100 from a refractory exponential distribution
set.seed(1102006,"Mersenne-Twister")
sampleSize <- 100
rate.true <- 20
rp.true <- 0.01
sampRE <- rrexp(sampleSize,rate=rate.true,rp=rp.true)
sampREmleRE <- rexpMLE(sampRE)
rbind(est = sampREmleRE$estimate,se = sampREmleRE$se,true = c(rate.true,rp.true))

## make a parametric boostrap to check the distribution of the deviance
nbReplicate <- 10000
system.time(
            devianceRE100 <- replicate(nbReplicate,{
              sampRE <- rrexp(sampleSize,rate=rate.true,rp=rp.true)
              sampREmleRE <- rexpMLE(sampRE)
              -2*sampREmleRE$r(rate.true,rp.true)
            }
                                       )
            )[3]

## Get 95 and 99\% confidence intervals for the QQ plot
ci <- sapply(1:nbReplicate,
                 function(idx) qchisq(qbeta(c(0.005,0.025,0.975,0.995),
                                            idx,
                                            nbReplicate-idx+1),
                                      df=2)
             )
## make QQ plot
X <- qchisq(ppoints(nbReplicate),df=2)
Y <- sort(devianceRE100)
X11()
plot(X,Y,type="n",
     xlab=expression(paste(chi[2]^2," quantiles")),
     ylab="MC quantiles",
     main="Deviance with true parameters after ML fit of refractory Poisson data",
     sub=paste("sample size:", sampleSize,"MC replicates:", nbReplicate)
     )
abline(a=0,b=1)
lines(X,ci[1,],lty=2)
lines(X,ci[2,],lty=2)
lines(X,ci[3,],lty=2)
lines(X,ci[4,],lty=2)
lines(X,Y,col=2)
}
}
\keyword{distribution}
\keyword{ts}
