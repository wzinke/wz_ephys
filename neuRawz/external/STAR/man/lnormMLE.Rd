\name{lnormMLE}
\alias{lnormMLE}
\title{Maximum Likelihood Parameter Estimation of a Log Normal Model with Possibly
  Censored Data}
\description{
  Estimate log normal model parameters by the maximum likelihood
  method using possibly censored data.
}
\usage{
lnormMLE(yi, ni = numeric(length(yi)) + 1,
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
  In the absence of censored data the ML estimates are available in
  closed form together with the Hessian matrix at
  the MLE. In presence of censored data an initial guess for the
  parameters is obtained using the uncensored data before maximizing the
  likelihood function to the full data set using \code{\link{optim}}
  with the \code{BFGS} method.

  In order to ensure good behavior of the numerical optimization
  routines, optimization is performed on the log of parameter
  \code{sdlog}.

  Standard errors are obtained from the inverse of the observed
  information matrix at the MLE. They are transformed to go from the log
  scale used by the optimization routine, when the latter is used (ie,
  for censored data) to the parameterization requested. 
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
}
\author{Christophe Pouzat  \email{christophe.pouzat@univ-paris5.fr} }
\note{
  The returned standard errors (component \code{se}) are valid in the asymptotic limit. You
  should plot contours using function \code{r} in the returned list and
  check that the contours are reasonably close to ellipses.
}
\seealso{\code{\link{Lognormal}},\code{\link{invgaussMLE}}}
\examples{
## Simulate sample of size 100 from a log normal
## distribution
set.seed(1102006,"Mersenne-Twister")
sampleSize <- 100
meanlog.true <- -2.4
sdlog.true <- 0.4
sampLN <- rlnorm(sampleSize,meanlog.true,sdlog.true)
sampLNmleLN <- lnormMLE(sampLN)
rbind(est = sampLNmleLN$estimate,se = sampLNmleLN$se,true = c(meanlog.true,sdlog.true))
## In the absence of censoring the MLE of the log normal is available in a
## closed form together with its variance (ie, the observed information matrix)
## we can check that we did not screw up at that stage by comparing the observed
## information matrix obtained numerically with the analytical one. To do that we
## use the MINUS log likelihood function returned by lnormMLE to get a numerical
## estimate
detailedFit <- optim(fn=sampLNmleLN$mll,
                     par=as.vector(c(sampLNmleLN$estimate[1],log(sampLNmleLN$estimate[2]))),
                     method="BFGS",
                     hessian=TRUE)
## We should not forget that the "mll" function uses the log of the sdlog parameter while
## the "se" component of sampLNmleLN list is expressed on the linear scale we must therefore
## transform one into the other as follows (Kalbfleisch, 1985, p71):
## if x = u and y = exp(v) and if we have the information matrix in term of
## u and v (that's the hessian component of list detailedFit above), we create matrix:
##      du/dx du/dy
## Q =
##      dv/dx dv/dy
## and we get I in term of x and y by the following matrix product:
## I(x,y) <- t(Q) \%*\% I(u,v) \%*\% Q
## In the present case:
##  du/dx = 1, du/dy = 0, dv/dx = 0, dv/dy = 1/exp(v)
## Therefore:
Q <- diag(c(1,1/exp(detailedFit$par[2])))
numericalI <- t(Q) \%*\% detailedFit$hessian \%*\% Q
seComp <- rbind(sampLNmleLN$se, sqrt(diag(solve(numericalI))))
colnames(seComp) <- c("meanlog","sdlog")
rownames(seComp) <- c("analytical", "numerical")
seComp
## We can check the relative differences between the 2
apply(seComp,2,function(x) abs(diff(x))/x[1])

\dontrun{
## Estimate the log relative likelihood on a grid to plot contours
MeanLog <- seq(sampLNmleLN$estimate[1]-4*sampLNmleLN$se[1],
               sampLNmleLN$estimate[1]+4*sampLNmleLN$se[1],
               sampLNmleLN$se[1]/10)
SdLog <- seq(sampLNmleLN$estimate[2]-4*sampLNmleLN$se[2],
             sampLNmleLN$estimate[2]+4*sampLNmleLN$se[2],
             sampLNmleLN$se[2]/10)
sampLNmleLNcontour <- sapply(MeanLog, function(mu) sapply(SdLog, function(s) sampLNmleLN$r(mu,s)))
## plot contours using a linear scale for the parameters
## draw four contours corresponding to the following likelihood ratios:
##  0.5, 0.1, Chi2 with 2 df and p values of 0.95 and 0.99
X11(width=12,height=6)
layout(matrix(1:2,ncol=2))
contour(MeanLog,SdLog,t(sampLNmleLNcontour),
        levels=c(log(c(0.5,0.1)),-0.5*qchisq(c(0.95,0.99),df=2)),
        labels=c("log(0.5)",
          "log(0.1)",
          "-1/2*P(Chi2=0.95)",
          "-1/2*P(Chi2=0.99)"),
        xlab=expression(mu),ylab=expression(sigma),
        main="Log Relative Likelihood Contours"
        )
points(sampLNmleLN$estimate[1],sampLNmleLN$estimate[2],pch=3)
points(meanlog.true,sdlog.true,pch=16,col=2)
## The contours are not really symmetrical about the MLE we can try to
## replot them using a log scale for the parameters to see if that improves
## the situation
contour(MeanLog,log(SdLog),t(sampLNmleLNcontour),
        levels=c(log(c(0.5,0.1)),-0.5*qchisq(c(0.95,0.99),df=2)),
        labels="",
        xlab=expression(mu),ylab=expression(log(sigma)),
        main="Log Relative Likelihood Contours",
        sub=expression(paste("log scale for parameter: ",sigma)))
points(sampLNmleLN$estimate[1],log(sampLNmleLN$estimate[2]),pch=3)
points(meanlog.true,log(sdlog.true),pch=16,col=2)

## make a parametric boostrap to check the distribution of the deviance
nbReplicate <- 10000
sampleSize <- 100
system.time(
            devianceLN100 <- replicate(nbReplicate,{
                             sampLN <- rlnorm(sampleSize,meanlog=meanlog.true,sdlog=sdlog.true)
                             sampLNmleLN <- lnormMLE(sampLN)
                             -2*sampLNmleLN$r(meanlog.true,sdlog.true)
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
Y <- sort(devianceLN100)
X11()
plot(X,Y,type="n",
     xlab=expression(paste(chi[2]^2," quantiles")),
     ylab="MC quantiles",
     main="Deviance with true parameters after ML fit of logNorm data",
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
