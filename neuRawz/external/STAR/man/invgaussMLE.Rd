\name{invgaussMLE}
\alias{invgaussMLE}
\title{Maximum Likelihood Parameter Estimation of an Inverse Gaussian Model with Possibly
  Censored Data}
\description{
  Estimate inverse Gaussian model parameters by the maximum likelihood
  method using possibly censored data. Two different parameterizations
  of the inverse Gaussian distribution can be used.
}
\usage{
invgaussMLE(yi, ni = numeric(length(yi)) + 1,
            si = numeric(length(yi)) + 1,
            parameterization = "sigma2")
}
\arguments{
  \item{yi}{vector of (possibly binned) observations or a
    \code{spikeTrain} object.}
  \item{ni}{vector of counts for each value of \code{yi}; default: \code{numeric(length(yi))+1}.}
  \item{si}{vector of counts of \emph{uncensored} observations for each
    value of \code{yi}; default: \code{numeric(length(yi))+1}.}
  \item{parameterization}{parameterization used, \code{"sigma2"}
    (default) of \code{"boundary"}.}
}
\details{
  The 2 different parameterizations of the inverse Gaussian distribution
  are discussed in the manual of \code{\link{dinvgauss}}.

  In the absence of censored data the ML estimates are available in
  closed form (Lindsey, 2004, p 212) together with the Hessian matrix at
  the MLE. In presence of censored data an initial guess for the
  parameters is obtained using the uncensored data before maximizing the
  likelihood function to the full data set using \code{\link{optim}}
  with the \code{BFGS} method. ML
  estimation is always performed with the \code{"sigma2"}
  parameterization. Parameters and variance-covariance matrix are
  transformed at the end if the \code{"boundary"} parameterization is
  requested.

  In order to ensure good behavior of the numerical optimization
  routines, optimization is performed on the log of the parameters
  (\code{mu} and \code{sigma2}).

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
    function using the log of the parameters.}
  \item{call}{the matched call.}
}
\references{
  Lindsey, J.K. (2004) \emph{Introduction to Applied Statistics: A
    Modelling Approach}. OUP.
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com} }
\note{
  The returned standard errors (component \code{se}) are valid in the asymptotic limit. You
  should plot contours using function \code{r} in the returned list and
  check that the contours are reasonably close to ellipses.
}
\seealso{\code{\link{dinvgauss}},\code{\link{lnormMLE}},\code{\link{gammaMLE}},\code{\link{weibullMLE}},\code{\link{llogisMLE}},\code{\link{rexpMLE}}.}
\examples{
## Simulate sample of size 100 from an inverse Gaussian
## distribution
set.seed(1102006,"Mersenne-Twister")
sampleSize <- 100
mu.true <- 0.075 
sigma2.true <- 3
sampleSize <- 100
sampIG <- rinvgauss(sampleSize,mu=mu.true,sigma2=sigma2.true)
## Make a maximum likelihood fit
sampIGmleIG <- invgaussMLE(sampIG)
## Compare estimates with actual values
rbind(est = coef(sampIGmleIG),se = sampIGmleIG$se,true = c(mu.true,sigma2.true))
## In the absence of censoring the MLE of the inverse Gaussian is available in a
## closed form together with its variance (ie, the observed information matrix)
## we can check that we did not screw up at that stage by comparing the observed
## information matrix obtained numerically with the analytical one. To do that we
## use the MINUS log likelihood function returned by invgaussMLE to get a numerical
## estimate
detailedFit <- optim(par=as.vector(log(sampIGmleIG$estimate)),
                     fn=sampIGmleIG$mll,
                     method="BFGS",
                     hessian=TRUE)
## We should not forget that the "mll" function uses the log of the parameters while
## the "se" component of sampIGmleIG list is expressed on the linear scale we must therefore
## transform one into the other as follows (Kalbfleisch, 1985, p71):
## if x = exp(u) and y = exp(v) and if we have the information matrix in term of
## u and v (that's the hessian component of list detailedFit above), we create matrix:
##      du/dx du/dy
## Q =
##      dv/dx dv/dy
## and we get I in term of x and y by the following matrix product:
## I(x,y) <- t(Q) \%*\% I(u,v) \%*\% Q
## In the present case:
##  du/dx = 1/exp(u), du/dy = 0, dv/dx = 0, dv/dy = 1/exp(v)
## Therefore:
Q <- diag(1/exp(detailedFit$par))
numericalI <- t(Q) \%*\% detailedFit$hessian \%*\% Q
seComp <- rbind(sampIGmleIG$se, sqrt(diag(solve(numericalI))))
colnames(seComp) <- c("mu","sigma2")
rownames(seComp) <- c("analytical", "numerical")
seComp
## We can check the relative differences between the 2
apply(seComp,2,function(x) abs(diff(x))/x[1])

\dontrun{
## Estimate the log relative likelihood on a grid to plot contours
Mu <- seq(coef(sampIGmleIG)[1]-4*sampIGmleIG$se[1],
          coef(sampIGmleIG)[1]+4*sampIGmleIG$se[1],
          sampIGmleIG$se[1]/10)
Sigma2 <- seq(coef(sampIGmleIG)[2]-4*sampIGmleIG$se[2],
              coef(sampIGmleIG)[2]+4*sampIGmleIG$se[2],
              sampIGmleIG$se[2]/10)
sampIGmleIGcontour <- sapply(Mu, function(mu) sapply(Sigma2, function(s2) sampIGmleIG$r(mu,s2)))
## plot contours using a linear scale for the parameters
## draw four contours corresponding to the following likelihood ratios:
##  0.5, 0.1, Chi2 with 2 df and p values of 0.95 and 0.99
X11(width=12,height=6)
layout(matrix(1:2,ncol=2))
contour(Mu,Sigma2,t(sampIGmleIGcontour),
        levels=c(log(c(0.5,0.1)),-0.5*qchisq(c(0.95,0.99),df=2)),
        labels=c("log(0.5)",
          "log(0.1)",
          "-1/2*P(Chi2=0.95)",
          "-1/2*P(Chi2=0.99)"),
        xlab=expression(mu),ylab=expression(sigma^2),
        main="Log Relative Likelihood Contours"
        )
points(coef(sampIGmleIG)[1],coef(sampIGmleIG)[2],pch=3)
points(mu.true,sigma2.true,pch=16,col=2)
## The contours are not really symmetrical about the MLE we can try to
## replot them using a log scale for the parameters to see if that improves
## the situation
contour(log(Mu),log(Sigma2),t(sampIGmleIGcontour),
        levels=c(log(c(0.5,0.1)),-0.5*qchisq(c(0.95,0.99),df=2)),
        labels="",
        xlab=expression(log(mu)),ylab=expression(log(sigma^2)),
        main="Log Relative Likelihood Contours",
        sub="log scale for the parameters")
points(log(coef(sampIGmleIG)[1]),log(coef(sampIGmleIG)[2]),pch=3)
points(log(mu.true),log(sigma2.true),pch=16,col=2)

## Even with the log scale the contours are not ellipsoidal, so let us compute profiles
## For that we are going to use the returned MINUS log likelihood function
logMuProfFct <- function(logMu,...) {
  myOpt <- optimise(function(x) sampIGmleIG$mll(c(logMu,x))+logLik(sampIGmleIG),...)
  as.vector(unlist(myOpt[c("objective","minimum")]))
}
logMuProfCI <- function(logMu,
                        CI,
                        a=logS2Seq[1],
                        b=logS2Seq[length(logS2Seq)]) logMuProfFct(logMu,c(a,b))[1] - qchisq(CI,1)/2

logS2ProfFct <- function(logS2,...) {
  myOpt <- optimise(function(x) sampIGmleIG$mll(c(x,logS2))+logLik(sampIGmleIG),...)
  as.vector(unlist(myOpt[c("objective","minimum")]))
}
logS2ProfCI <- function(logS2, CI,
                        a=logMuSeq[1],
                        b=logMuSeq[length(logMuSeq)]) logS2ProfFct(logS2,c(a,b))[1] - qchisq(CI,1)/2


## We compute profiles (on the log scale) eploxing +/- 3 times
## the se about the MLE
logMuSE <- sqrt(diag(solve(detailedFit$hessian)))[1]
logMuSeq <- seq(log(coef(sampIGmleIG)[1])-3*logMuSE,
                log(coef(sampIGmleIG)[1])+3*logMuSE,
                logMuSE/10)
logS2SE <- sqrt(diag(solve(detailedFit$hessian)))[2]
logS2Seq <- seq(log(coef(sampIGmleIG)[2])-3*logS2SE,
                log(coef(sampIGmleIG)[2])+3*logS2SE,
                logS2SE/10)
logMuProf <- sapply(logMuSeq,logMuProfFct,
                    lower=logS2Seq[1],
                    upper=logS2Seq[length(logS2Seq)])
## Get 95% and 99% CI
logMuCI95 <- c(uniroot(logMuProfCI,
                       interval=c(logMuSeq[1],log(coef(sampIGmleIG)[1])),
                       CI=0.95)$root,
               uniroot(logMuProfCI,
                       interval=c(log(coef(sampIGmleIG)[1]),logMuSeq[length(logMuSeq)]),
                       CI=0.95)$root
               )
logMuCI99 <- c(uniroot(logMuProfCI,
                       interval=c(logMuSeq[1],log(coef(sampIGmleIG)[1])),
                       CI=0.99)$root,
               uniroot(logMuProfCI,
                       interval=c(log(coef(sampIGmleIG)[1]),logMuSeq[length(logMuSeq)]),
                       CI=0.99)$root
               )

logS2Prof <- sapply(logS2Seq,logS2ProfFct,
                    lower=logMuSeq[1],
                    upper=logMuSeq[length(logMuSeq)])
## Get 95% and 99% CI
logS2CI95 <- c(uniroot(logS2ProfCI,
                       interval=c(logS2Seq[1],log(coef(sampIGmleIG)[2])),
                       CI=0.95)$root,
               uniroot(logS2ProfCI,
                       interval=c(log(coef(sampIGmleIG)[2]),logS2Seq[length(logS2Seq)]),
                       CI=0.95)$root
               )
logS2CI99 <- c(uniroot(logS2ProfCI,
                       interval=c(logS2Seq[1],log(coef(sampIGmleIG)[2])),
                       CI=0.99)$root,
               uniroot(logS2ProfCI,
                       interval=c(log(coef(sampIGmleIG)[2]),logS2Seq[length(logS2Seq)]),
                       CI=0.99)$root
               )


## Add profiles to the previous plot
lines(logMuSeq,logMuProf[2,],col=2,lty=2)
lines(logS2Prof[2,],logS2Seq,col=2,lty=2)

## We can now check the deviations of the (profiled) deviances
## from the asymptotic parabolic curves
X11()
layout(matrix(1:4,nrow=2))
oldpar <- par(mar=c(4,4,2,1))
logMuSeqOffset <- logMuSeq-log(coef(sampIGmleIG)[1])
logMuVar <- diag(solve(detailedFit$hessian))[1]
plot(logMuSeq,2*logMuProf[1,],type="l",xlab=expression(log(mu)),ylab="Deviance")
lines(logMuSeq,logMuSeqOffset^2/logMuVar,col=2)
points(log(coef(sampIGmleIG)[1]),0,pch=3)
abline(h=0)
abline(h=qchisq(0.95,1),lty=2)
abline(h=qchisq(0.99,1),lty=2)
lines(rep(logMuCI95[1],2),c(0,qchisq(0.95,1)),lty=2)
lines(rep(logMuCI95[2],2),c(0,qchisq(0.95,1)),lty=2)
lines(rep(logMuCI99[1],2),c(0,qchisq(0.99,1)),lty=2)
lines(rep(logMuCI99[2],2),c(0,qchisq(0.99,1)),lty=2)
## We can also "linearize" this last graph
plot(logMuSeq,
     sqrt(2*logMuProf[1,])*sign(logMuSeqOffset),
     type="l",
     xlab=expression(log(mu)),
     ylab=expression(paste("signed ",sqrt(Deviance)))
     )
lines(logMuSeq,
      sqrt(logMuSeqOffset^2/logMuVar)*sign(logMuSeqOffset),
      col=2)
points(log(coef(sampIGmleIG)[1]),0,pch=3)

logS2SeqOffset <- logS2Seq-log(coef(sampIGmleIG)[2])
logS2Var <- diag(solve(detailedFit$hessian))[2]
plot(logS2Seq,2*logS2Prof[1,],type="l",xlab=expression(log(sigma^2)),ylab="Deviance")
lines(logS2Seq,logS2SeqOffset^2/logS2Var,col=2)
points(log(coef(sampIGmleIG)[2]),0,pch=3)
abline(h=0)
abline(h=qchisq(0.95,1),lty=2)
abline(h=qchisq(0.99,1),lty=2)
lines(rep(logS2CI95[1],2),c(0,qchisq(0.95,1)),lty=2)
lines(rep(logS2CI95[2],2),c(0,qchisq(0.95,1)),lty=2)
lines(rep(logS2CI99[1],2),c(0,qchisq(0.99,1)),lty=2)
lines(rep(logS2CI99[2],2),c(0,qchisq(0.99,1)),lty=2)
## We can also "linearize" this last graph
plot(logS2Seq,
     sqrt(2*logS2Prof[1,])*sign(logS2SeqOffset),
     type="l",
     xlab=expression(log(sigma^2)),
     ylab=expression(paste("signed ",sqrt(Deviance)))
     )
lines(logS2Seq,
      sqrt(logS2SeqOffset^2/logS2Var)*sign(logS2SeqOffset),
      col=2)
points(log(coef(sampIGmleIG)[2]),0,pch=3)
par(oldpar)

## make a parametric boostrap to check the distribution of the deviance
nbReplicate <- 1000 #10000
sampleSize <- 100
system.time(
devianceIG100 <- lapply(1:nbReplicate,
                        function(idx) {
                          if ((idx %% 10 - 1) == 0) cat(paste("Doing now iteration:",idx,"\n"))
                          sampIG <- rinvgauss(sampleSize,mu=mu.true,sigma2=sigma2.true)
                          sampIGmleIG <- invgaussMLE(sampIG)
                          Deviance <- -2*sampIGmleIG$r(mu.true,sigma2.true)
                          logPara <- log(coef(sampIGmleIG))
                          logParaSE <- sampIGmleIG$se/coef(sampIGmleIG)
                          intervalMu <- function(n) c(-n,n)*logParaSE[1]+logPara[1]
                          intervalS2 <- function(n) c(-n,n)*logParaSE[2]+logPara[2]
                          logMuProfFct <- function(logMu,...) {
                            optimise(function(x)
                                     sampIGmleIG$mll(c(logMu,x))+logLik(sampIGmleIG),...)$objective
                          }
                          logMuProfCI <- function(logMu,
                                                  CI,
                                                  a=intervalS2(4)[1],
                                                  b=intervalS2(4)[2])
                            logMuProfFct(logMu,c(a,b)) - qchisq(CI,1)/2
                          
                          logS2ProfFct <- function(logS2,...) {
                            optimise(function(x)
                                     sampIGmleIG$mll(c(x,logS2))+logLik(sampIGmleIG),...)$objective
                          }
                          logS2ProfCI <- function(logS2, CI,
                                                  a=intervalMu(4)[1],
                                                  b=intervalMu(4)[2])
                            logS2ProfFct(logS2,c(a,b)) - qchisq(CI,1)/2
                          
                          factor <- 4
                          while((logMuProfCI(intervalMu(factor)[2],0.99) *
                                 logMuProfCI(logPara[1],0.99) >= 0) ||
                                (logMuProfCI(intervalMu(factor)[1],0.99) *
                                 logMuProfCI(logPara[1],0.99) >= 0)
                                ) factor <- factor+1
                          ##browser()
                          logMuCI95 <- c(uniroot(logMuProfCI,
                                                 interval=c(intervalMu(factor)[1],logPara[1]),
                                                 CI=0.95)$root,
                                         uniroot(logMuProfCI,
                                                 interval=c(logPara[1],intervalMu(factor)[2]),
                                                 CI=0.95)$root
                                         )
                          logMuCI99 <- c(uniroot(logMuProfCI,
                                                 interval=c(intervalMu(factor)[1],logPara[1]),
                                                 CI=0.99)$root,
                                         uniroot(logMuProfCI,
                                                 interval=c(logPara[1],intervalMu(factor)[2]),
                                                 CI=0.99)$root
                                         )
                          factor <- 4
                          while((logS2ProfCI(intervalS2(factor)[2],0.99) *
                                 logS2ProfCI(logPara[2],0.99) >= 0) ||
                                (logS2ProfCI(intervalS2(factor)[1],0.99) *
                                 logS2ProfCI(logPara[2],0.99) >= 0)
                                ) factor <- factor+1
                          logS2CI95 <- c(uniroot(logS2ProfCI,
                                                 interval=c(intervalS2(factor)[1],logPara[2]),
                                                 CI=0.95)$root,
                                         uniroot(logS2ProfCI,
                                                    interval=c(logPara[2],intervalS2(factor)[2]),
                                                 CI=0.95)$root
                                         )
                          logS2CI99 <- c(uniroot(logS2ProfCI,
                                                 interval=c(intervalS2(factor)[1],logPara[2]),
                                                 CI=0.99)$root,
                                         uniroot(logS2ProfCI,
                                                 interval=c(logPara[2],intervalS2(factor)[2]),
                                                 CI=0.99)$root
                                         )
                          list(deviance=Deviance,
                               logMuCI95=logMuCI95,
                               logMuNorm95=qnorm(c(0.025,0.975),logPara[1],logParaSE[1]),
                               logMuCI99=logMuCI99,
                               logMuNorm99=qnorm(c(0.005,0.995),logPara[1],logParaSE[1]),
                               logS2CI95=logS2CI95,
                               logS2Norm95=qnorm(c(0.025,0.975),logPara[2],logParaSE[2]),
                               logS2CI99=logS2CI99,
                               logS2Norm99=qnorm(c(0.005,0.995),logPara[2],logParaSE[2]))
                        }
                        )
            )[3]
## Find out how many times the true parameters was within the computed CIs
nLogMuCI95 <- sum(sapply(devianceIG100,
                         function(l) l$logMuCI95[1] <= log(mu.true)  &&
                         log(mu.true)<= l$logMuCI95[2]
                         )
                  )
nLogMuNorm95 <- sum(sapply(devianceIG100,
                           function(l) l$logMuNorm95[1] <= log(mu.true)  &&
                           log(mu.true)<= l$logMuNorm95[2]
                           )
                    )
nLogMuCI99 <- sum(sapply(devianceIG100,
                         function(l) l$logMuCI99[1] <= log(mu.true)  &&
                         log(mu.true)<= l$logMuCI99[2]
                         )
                  )
nLogMuNorm99 <- sum(sapply(devianceIG100,
                           function(l) l$logMuNorm99[1] <= log(mu.true)  &&
                           log(mu.true)<= l$logMuNorm99[2]
                           )
                    )
## Check if these counts are compatible with the nominal CIs
c(prof95Mu=nLogMuCI95,norm95Mu=nLogMuNorm95)
qbinom(c(0.005,0.995),nbReplicate,0.95)
c(prof95Mu=nLogMuCI99,norm95Mu=nLogMuNorm99)
qbinom(c(0.005,0.995),nbReplicate,0.99)

nLogS2CI95 <- sum(sapply(devianceIG100,
                         function(l) l$logS2CI95[1] <= log(sigma2.true)  &&
                         log(sigma2.true)<= l$logS2CI95[2]
                         )
                  )
nLogS2Norm95 <- sum(sapply(devianceIG100,
                           function(l) l$logS2Norm95[1] <= log(sigma2.true)  &&
                           log(sigma2.true)<= l$logS2Norm95[2]
                           )
                    )
nLogS2CI99 <- sum(sapply(devianceIG100,
                         function(l) l$logS2CI99[1] <= log(sigma2.true)  &&
                         log(sigma2.true)<= l$logS2CI99[2]
                         )
                  )
nLogS2Norm99 <- sum(sapply(devianceIG100,
                           function(l) l$logS2Norm99[1] <= log(sigma2.true)  &&
                           log(sigma2.true)<= l$logS2Norm99[2]
                           )
                    )
## Check if these counts are compatible with the nominal CIs
c(prof95S2=nLogS2CI95,norm95S2=nLogS2Norm95)
qbinom(c(0.005,0.995),nbReplicate,0.95)
c(prof95S2=nLogS2CI99,norm95S2=nLogS2Norm99)
qbinom(c(0.005,0.995),nbReplicate,0.99)


## Get 95 and 99\% confidence intervals for the QQ plot
ci <- sapply(1:nbReplicate,
                 function(idx) qchisq(qbeta(c(0.005,0.025,0.975,0.995),
                                            idx,
                                            nbReplicate-idx+1),
                                      df=2)
             )
## make QQ plot
X <- qchisq(ppoints(nbReplicate),df=2)
Y <- sort(sapply(devianceIG100,function(l) l$deviance))
X11()
plot(X,Y,type="n",
     xlab=expression(paste(chi[2]^2," quantiles")),
     ylab="MC quantiles",
     main="Deviance with true parameters after ML fit of IG data",
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
