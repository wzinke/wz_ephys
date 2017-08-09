\name{weibullMLE}
\alias{weibullMLE}
\title{Maximum Likelihood Parameter Estimation of a Weibull Model with Possibly
  Censored Data}
\description{
  Estimate Weibull model parameters by the maximum likelihood
  method using possibly censored data.
}
\usage{
weibullMLE(yi, ni = numeric(length(yi)) + 1,
           si = numeric(length(yi)) + 1, shape.min = 0.05, shape.max = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yi}{vector of (possibly binned) observations or a
    \code{spikeTrain} object.}
  \item{ni}{vector of counts for each value of \code{yi}; default: \code{numeric(length(yi))+1}.}
  \item{si}{vector of counts of \emph{uncensored} observations for each
    value of \code{yi}; default: \code{numeric(length(yi))+1}.}
  \item{shape.min}{numeric, the inital guess of the minimal possible
    value of the \code{shape} parameter, used by \code{optimise}.}
  \item{shape.max}{numeric, the inital guess of the maximal possible
    value of the \code{shape} parameter, used by \code{optimise}.}
}
\details{
  There is no closed form expression for the MLE of a Weibull distribution. The
  numerical method implemented here uses the profile likelihood
  described by Kalbfleisch (1985) pp 56-58.

  In order to ensure good behavior of the numerical optimization
  routines, optimization is performed on the log of the parameters
  (\code{shape} and \code{scale}).

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
  Kalbfleisch, J. G. (1985) \emph{Probability and Statistical
    Inference. Volume 2: Statistical Inference}. Springer-Verlag.

  Lindsey, J.K. (2004) \emph{Introduction to Applied Statistics: A
    Modelling Approach}. OUP.
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com} }
\note{
  The returned standard errors (component \code{se}) are valid in the asymptotic limit. You
  should plot contours using function \code{r} in the returned list and
  check that the contours are reasonably close to ellipses.
}
\seealso{\code{\link{Weibull}},
  \code{\link{invgaussMLE}},
  \code{\link{lnormMLE}},
  \code{\link{gammaMLE}}}
\examples{
\dontrun{
## Simulate sample of size 100 from a weibull distribution
set.seed(1102006,"Mersenne-Twister")
sampleSize <- 100
shape.true <- 2.5
scale.true <- 0.085
sampWB <- rweibull(sampleSize,shape=shape.true,scale=scale.true)
sampWBmleWB <- weibullMLE(sampWB)
rbind(est = sampWBmleWB$estimate,se = sampWBmleWB$se,true = c(shape.true,scale.true))

## Estimate the log relative likelihood on a grid to plot contours
Shape <- seq(sampWBmleWB$estimate[1]-4*sampWBmleWB$se[1],
               sampWBmleWB$estimate[1]+4*sampWBmleWB$se[1],
               sampWBmleWB$se[1]/10)
Scale <- seq(sampWBmleWB$estimate[2]-4*sampWBmleWB$se[2],
             sampWBmleWB$estimate[2]+4*sampWBmleWB$se[2],
             sampWBmleWB$se[2]/10)
sampWBmleWBcontour <- sapply(Shape, function(sh) sapply(Scale, function(sc) sampWBmleWB$r(sh,sc)))
## plot contours using a linear scale for the parameters
## draw four contours corresponding to the following likelihood ratios:
##  0.5, 0.1, Chi2 with 2 df and p values of 0.95 and 0.99
X11(width=12,height=6)
layout(matrix(1:2,ncol=2))
contour(Shape,Scale,t(sampWBmleWBcontour),
        levels=c(log(c(0.5,0.1)),-0.5*qchisq(c(0.95,0.99),df=2)),
        labels=c("log(0.5)",
          "log(0.1)",
          "-1/2*P(Chi2=0.95)",
          "-1/2*P(Chi2=0.99)"),
        xlab="shape",ylab="scale",
        main="Log Relative Likelihood Contours"
        )
points(sampWBmleWB$estimate[1],sampWBmleWB$estimate[2],pch=3)
points(shape.true,scale.true,pch=16,col=2)
## The contours are not really symmetrical about the MLE we can try to
## replot them using a log scale for the parameters to see if that improves
## the situation
contour(log(Shape),log(Scale),t(sampWBmleWBcontour),
        levels=c(log(c(0.5,0.1)),-0.5*qchisq(c(0.95,0.99),df=2)),
        labels="",
        xlab="log(shape)",ylab="log(scale)",
        main="Log Relative Likelihood Contours",
        sub="log scale for the parameters")
points(log(sampWBmleWB$estimate[1]),log(sampWBmleWB$estimate[2]),pch=3)
points(log(shape.true),log(scale.true),pch=16,col=2)

## make a parametric boostrap to check the distribution of the deviance
nbReplicate <- 10000
sampleSize <- 100
system.time(
            devianceWB100 <- replicate(nbReplicate,{
              sampWB <- rweibull(sampleSize,shape=shape.true,scale=scale.true)
              sampWBmleWB <- weibullMLE(sampWB)
              -2*sampWBmleWB$r(shape.true,scale.true)
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
Y <- sort(devianceWB100)
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
