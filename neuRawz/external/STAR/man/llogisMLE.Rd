\name{llogisMLE}
\alias{llogisMLE}
\title{Maximum Likelihood Parameter Estimation of a Log Logistic Model with Possibly
  Censored Data}
\description{
  Estimate log logistic model parameters by the maximum likelihood
  method using possibly censored data.
}
\usage{
llogisMLE(yi, ni = numeric(length(yi)) + 1,
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
  The MLE for the log logistic is not available in closed formed and
  is therefore obtained numerically obtained by calling
  \code{\link{optim}} with the \code{BFGS} method.

  In order to ensure good behavior of the numerical optimization
  routines, optimization is performed on the log of parameter
  \code{scale}.

  Standard errors are obtained from the inverse of the observed
  information matrix at the MLE. They are transformed to go from the log
  scale used by the optimization routine to the requested parameterization. 
}
\value{
  A list of class \code{durationFit} with the following components:
  \item{estimate}{the estimated parameters, a named vector.}
  \item{se}{the standard errors, a named vector.}
  \item{logLik}{the log likelihood at maximum.}
  \item{r}{a function returning the log of the relative likelihood function.}
  \item{mll}{a function returning the opposite of the log likelihood
    function using the log of parameter \code{sdlog}.}
  \item{call}{the matched call.}
}
\references{
  Lindsey, J.K. (2004) \emph{Introduction to Applied Statistics: A
    Modelling Approach}. OUP.

  Lindsey, J.K. (2004) \emph{The Statistical Analysis of Stochastic
    Processes in Time}. CUP.
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com} }
\note{
  The returned standard errors (component \code{se}) are valid in the asymptotic limit. You
  should plot contours using function \code{r} in the returned list and
  check that the contours are reasonably close to ellipses.
}
\seealso{
  \code{\link{dllogis}},
  \code{\link{invgaussMLE}},
  \code{\link{gammaMLE}},
  \code{\link{weibullMLE}},
  \code{\link{rexpMLE}},
  \code{\link{lnormMLE}}
}
\examples{
\dontrun{
## Simulate sample of size 100 from a log logisitic
## distribution
set.seed(1102006,"Mersenne-Twister")
sampleSize <- 100
location.true <- -2.7
scale.true <- 0.025
sampLL <- rllogis(sampleSize,location=location.true,scale=scale.true)
sampLLmleLL <- llogisMLE(sampLL)
rbind(est = sampLLmleLL$estimate,se = sampLLmleLL$se,true = c(location.true,scale.true))

## Estimate the log relative likelihood on a grid to plot contours
Loc <- seq(sampLLmleLL$estimate[1]-4*sampLLmleLL$se[1],
               sampLLmleLL$estimate[1]+4*sampLLmleLL$se[1],
               sampLLmleLL$se[1]/10)
Scale <- seq(sampLLmleLL$estimate[2]-4*sampLLmleLL$se[2],
             sampLLmleLL$estimate[2]+4*sampLLmleLL$se[2],
             sampLLmleLL$se[2]/10)
sampLLmleLLcontour <- sapply(Loc, function(m) sapply(Scale, function(s) sampLLmleLL$r(m,s)))
## plot contours using a linear scale for the parameters
## draw four contours corresponding to the following likelihood ratios:
##  0.5, 0.1, Chi2 with 2 df and p values of 0.95 and 0.99
X11(width=12,height=6)
layout(matrix(1:2,ncol=2))
contour(Loc,Scale,t(sampLLmleLLcontour),
        levels=c(log(c(0.5,0.1)),-0.5*qchisq(c(0.95,0.99),df=2)),
        labels=c("log(0.5)",
          "log(0.1)",
          "-1/2*P(Chi2=0.95)",
          "-1/2*P(Chi2=0.99)"),
        xlab="Location",ylab="Scale",
        main="Log Relative Likelihood Contours"
        )
points(sampLLmleLL$estimate[1],sampLLmleLL$estimate[2],pch=3)
points(location.true,scale.true,pch=16,col=2)
## The contours are not really symmetrical about the MLE we can try to
## replot them using a log scale for the parameters to see if that improves
## the situation
contour(Loc,log(Scale),t(sampLLmleLLcontour),
        levels=c(log(c(0.5,0.1)),-0.5*qchisq(c(0.95,0.99),df=2)),
        labels="",
        xlab="log(Location)",ylab="log(Scale)",
        main="Log Relative Likelihood Contours",
        sub="log scale for parameter: scale")
points(sampLLmleLL$estimate[1],log(sampLLmleLL$estimate[2]),pch=3)
points(location.true,log(scale.true),pch=16,col=2)

## make a parametric boostrap to check the distribution of the deviance
nbReplicate <- 10000
sampleSize <- 100
system.time(
            devianceLL100 <- replicate(nbReplicate,{
              sampLL <- rllogis(sampleSize,location=location.true,scale=scale.true)
              sampLLmleLL <- llogisMLE(sampLL)
              -2*sampLLmleLL$r(location.true,scale.true)
            }
                                       )
            )[3]

## Get 95 and 99% confidence intervals for the QQ plot
ci <- sapply(1:nbReplicate,
                 function(idx) qchisq(qbeta(c(0.005,0.025,0.975,0.995),
                                            idx,
                                            nbReplicate-idx+1),
                                      df=2)
             )
## make QQ plot
X <- qchisq(ppoints(nbReplicate),df=2)
Y <- sort(devianceLL100)
X11()
plot(X,Y,type="n",
     xlab=expression(paste(chi[2]^2," quantiles")),
     ylab="MC quantiles",
     main="Deviance with true parameters after ML fit of log logistic data",
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
