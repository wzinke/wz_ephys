\name{gammaMLE}
\alias{gammaMLE}
\title{Maximum Likelihood Parameter Estimation of a Gamma Model with Possibly
  Censored Data}
\description{
  Estimate Gamma model parameters by the maximum likelihood
  method using possibly censored data. Two different parameterizations
  of the Gamma distribution can be used.
}
\usage{
gammaMLE(yi, ni = numeric(length(yi)) + 1,
         si = numeric(length(yi)) + 1, scale = TRUE)
}
\arguments{
  \item{yi}{vector of (possibly binned) observations or a
    \code{spikeTrain} object.}
  \item{ni}{vector of counts for each value of \code{yi}; default: \code{numeric(length(yi))+1}.}
  \item{si}{vector of counts of \emph{uncensored} observations for each
    value of \code{yi}; default: \code{numeric(length(yi))+1}.}
  \item{scale}{logical should the scale (\code{TRUE}) or the rate
    parameterization (\code{FALSE}) be used?}
}
\details{
  There is no closed form expression for the MLE of a Gamma distribution. The
  numerical method implemented here uses the profile likelihood
  described by Monahan (2001) pp 210-216.

  In order to ensure good behavior of the numerical optimization
  routines, optimization is performed on the log of the parameters
  (\code{shape} and \code{scale} or \code{rate}).

  Standard errors are obtained from the inverse of the observed
  information matrix at the MLE. They are transformed to go from the log
  scale used by the optimization routine to the parameterization requested.
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
\references{
  Monahan, J. F. (2001) \emph{Numerical Methods of Statistics}. CUP.
  
  Lindsey, J.K. (2004) \emph{Introduction to Applied Statistics: A
    Modelling Approach}. OUP.
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com} }
\note{
  The returned standard errors (component \code{se}) are valid in the asymptotic limit. You
  should plot contours using function \code{r} in the returned list and
  check that the contours are reasonably close to ellipses.
}
\seealso{
  \code{\link{GammaDist}},
  \code{\link{invgaussMLE}},
  \code{\link{lnormMLE}}}
\examples{
\dontrun{
## Simulate sample of size 100 from a gamma distribution
set.seed(1102006,"Mersenne-Twister")
sampleSize <- 100
shape.true <- 6
scale.true <- 0.012
sampGA <- rgamma(sampleSize,shape=shape.true,scale=scale.true)
sampGAmleGA <- gammaMLE(sampGA)
rbind(est = sampGAmleGA$estimate,se = sampGAmleGA$se,true = c(shape.true,scale.true))

## Estimate the log relative likelihood on a grid to plot contours
Shape <- seq(sampGAmleGA$estimate[1]-4*sampGAmleGA$se[1],
               sampGAmleGA$estimate[1]+4*sampGAmleGA$se[1],
               sampGAmleGA$se[1]/10)
Scale <- seq(sampGAmleGA$estimate[2]-4*sampGAmleGA$se[2],
             sampGAmleGA$estimate[2]+4*sampGAmleGA$se[2],
             sampGAmleGA$se[2]/10)
sampGAmleGAcontour <- sapply(Shape, function(sh) sapply(Scale, function(sc) sampGAmleGA$r(sh,sc)))
## plot contours using a linear scale for the parameters
## draw four contours corresponding to the following likelihood ratios:
##  0.5, 0.1, Chi2 with 2 df and p values of 0.95 and 0.99
X11(width=12,height=6)
layout(matrix(1:2,ncol=2))
contour(Shape,Scale,t(sampGAmleGAcontour),
        levels=c(log(c(0.5,0.1)),-0.5*qchisq(c(0.95,0.99),df=2)),
        labels=c("log(0.5)",
          "log(0.1)",
          "-1/2*P(Chi2=0.95)",
          "-1/2*P(Chi2=0.99)"),
        xlab="shape",ylab="scale",
        main="Log Relative Likelihood Contours"
        )
points(sampGAmleGA$estimate[1],sampGAmleGA$estimate[2],pch=3)
points(shape.true,scale.true,pch=16,col=2)
## The contours are not really symmetrical about the MLE we can try to
## replot them using a log scale for the parameters to see if that improves
## the situation
contour(log(Shape),log(Scale),t(sampGAmleGAcontour),
        levels=c(log(c(0.5,0.1)),-0.5*qchisq(c(0.95,0.99),df=2)),
        labels="",
        xlab="log(shape)",ylab="log(scale)",
        main="Log Relative Likelihood Contours",
        sub="log scale for the parameters")
points(log(sampGAmleGA$estimate[1]),log(sampGAmleGA$estimate[2]),pch=3)
points(log(shape.true),log(scale.true),pch=16,col=2)

## make a parametric boostrap to check the distribution of the deviance
nbReplicate <- 10000
sampleSize <- 100
system.time(
            devianceGA100 <- replicate(nbReplicate,{
                             sampGA <- rgamma(sampleSize,shape=shape.true,scale=scale.true)
                             sampGAmleGA <- gammaMLE(sampGA)
                             -2*sampGAmleGA$r(shape.true,scale.true)
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
Y <- sort(devianceGA100)
X11()
plot(X,Y,type="n",
     xlab=expression(paste(chi[2]^2," quantiles")),
     ylab="MC quantiles",
     main="Deviance with true parameters after ML fit of gamma data",
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
