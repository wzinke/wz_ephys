\name{mkGLMdf}
\alias{mkGLMdf}
\title{ Formats (lists of) spikeTrain and repeatedTrain Objects into
  Data Frame for use in glm, mgcv and gam }
\description{
  Given a \code{spikeTrain} or a \code{repeatedTrain} objects or a
  \code{list} of any of those two, \code{mkGLMdf} generates a
  \code{data.frame}, by discretizing time, allowing \code{glm}, \code{gss} and
  \code{gam} to be used with the \code{poisson} or \code{binomial}
  family to fit the spike trains.
}
\usage{
mkGLMdf(obj, delta, lwr, upr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{a \code{spikeTrain} or a \code{repeatedTrain} objects or a
  \code{list} of any of those two.}
  \item{delta}{the bin size used for time discretization (in s).}
  \item{lwr}{the time (in s) at which the recording window starts. If
    \code{missing} a value is obtained
    using the \code{\link{floor}} of the smallest spike time.}
  \item{upr}{the time (in s) at which the recording window ends. If
    \code{missing} a value is obtained
    using the \code{\link{ceiling}} of the largest spike time.}
  
}
\details{
  The construction of the returned list is very clearly explained in Jim
  Lindsey's paper (1995). The idea has been used several time in the
  field: Brillinger (1988), Kass and Ventura (2001), Truccolo et al (2005). 
}
\value{
  A \code{\link{data.frame}} with the following variables:
  \item{event }{an integer presence (1) or absence (0) of an event from
    a given neuron in the given bin.}
  \item{time }{time at bin center.}
  \item{neuron }{a factor giving the neuron to which this row of the
    data frame refers.}
  \item{lN.x }{a \code{numeric}. \code{x} takes value
    1, 2, ..., number of neurons present in \code{obj}. The time to the
    last event of the corresponding neuron.}

  The list has also few attributes: \code{lwr}, the start of the
  recording window; \code{upr}, the end of the recording window;
  \code{delta}, the bin width; \code{call}, the call used to generate
  the list.
}
\references{
  Lindsey, J. K. (1995) Fitting Parametric Counting Processes by Using
  Log-Linear Models \emph{Applied Statistics} \bold{44}: 201--212.

  Brillinger, D. R. (1988) Maximum likelihood analysis of spike trains of
  interacting nerve cells \emph{Biol Cybern} \bold{59}: 189--200.

  Kass, Robert E. and Ventura, Val\'erie (2001) A spike-train
  probability model \emph{Neural Comput.} \bold{13}: 1713--1720.

  Truccolo, W., Eden, U. T., Fellows, M. R., Donoghue, J. P. and Brown,
  E. N. (2005) A Point Process Framework for Relating Neural Spiking
  Activity to Spiking History, Neural Ensemble and Extrinsic Covariate
  Effects \emph{J Neurophysiol} \bold{93}:
  1074--1089. \url{http://jn.physiology.org/cgi/content/abstract/93/2/1074} 
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com}}
\note{
  See the example bellow to get an idea of what to do with the returned list.
}
\seealso{
  \code{\link{data.frame}},
  \code{\link{glm}},
  \code{\link[gss]{gssanova}},
  \code{\link{mgcv}},
  \code{\link{as.spikeTrain}},
  \code{\link{as.repeatedTrain}}
}
\examples{
\dontrun{
## Analysis of a "simple" spontaneous train
## load the data
data(e060824spont)
## create a data frame using the 1st neuron
DFA <- subset(mkGLMdf(e060824spont,0.004,0,59),neuron==1)
## Add the previous ISI to the data frame
DFA <- within(DFA,i1 <- isi(DFA,lag=1))
DFA <- DFA[complete.cases(DFA),]
## estimate the "map to uniform" functions for 2 variables:
## 1) the elapsed time since the last spike (lN.1)
## 2) the previous insterspike interval
## Do this estimation on the first half of the set
m2u1 <- mkM2U(DFA,"lN.1",0,29)
m2ui <- mkM2U(DFA,"i1",0,29,maxiter=200)
## create "mapped" variables
DFA <- within(DFA,e1t <- m2u1(lN.1))
DFA <- within(DFA,i1t <- m2ui(i1))
## split the data in 2 parts, one for the fit, the other for the test
DFAe <- subset(DFA, time <= 29)
DFAl <- subset(DFA, time > 29)
## fit an additive model with gssanova
m1.fit <- gssanova(event ~ e1t + i1t,
                   data = DFAe,
                   family="binomial",
                   seed=20061001)
## test the model by "time transforming" the late part
tt.l <- m1.fit \%tt\% DFAl
tt.l.summary <- summary(tt.l)
tt.l.summary
plot(tt.l.summary,which=c(1,2,4,6))

## Start with simulatd data #####
## Use thinning method and for that define a couple
## of functions

## expDecay gives an exponentially decaying
## synaptic effect followin a presynpatic spike.
## All the pre-synaptic spikes between "now" (argument
## t) and the previous spike of the post-synaptic
## neuron have an effect (and the summation is linear)
expDecay <- function(t,preT,last,
                     delay=0.002,tau=0.015) {
  
  if (missing(last)) good <- (preT+delay) < t
  else good <- last < preT & (preT+delay) < t
  if (sum(good) == 0) return(0)
  preS <- preT[good]
  preS <- t-preS-delay
  sum(exp(-preS/tau))

}

## Same as expDecay except that the effect is pusle like
pulseFF <- function(t,preT,last,
                    delay=0.005,duration=0.01) {
  if (missing(last)) good <- t-duration < (preT+delay) & (preT+delay) < t
  else good <- t-duration < (preT+delay) & last < preT & (preT+delay) < t
  sum(good)
}

## The work horse. Given a pre-synaptic train (preT),
## a duration, lognormal parameters and a presynaptic
## effect fucntion, mkPostTrain simulates a log-linear
## post-synaptic train using the thinning method
mkPostTrain <- function(preT,
                        duration=60,
                        meanlog=-2.4,
                        sdlog=0.4,
                        preFF=expDecay,
                        beta=log(5),
                        maxCI=30,
                        ...) {

  nuRest <- exp(-meanlog-0.5*sdlog^2)
  poissonRest <- nuRest*ifelse(beta>0,exp(beta),1)
  ciRest <- function(t) nuRest*exp(beta*preFF(t,preT,...))

  poissonNext <- maxCI*ifelse(beta>0,exp(beta),1)
  ci <- function(t,tLast) hlnorm(t-tLast,meanlog,sdlog)*exp(beta*preFF(t,preT,tLast,...))

  vLength <- poissonRest*300
  result <- numeric(vLength)
  currentTime <- 0
  lastTime <- 0
  eventIdx <- 1

  nextTime <- function(currentTime,lastTime) {
    if (currentTime > 0) {
      currentTime <- currentTime + rexp(1,poissonNext)
      ciRatio <- ci(currentTime,lastTime)/poissonNext
      if (ciRatio > 1) stop("Problem with thinning.")
      while (runif(1) > ciRatio) {
        currentTime <- currentTime + rexp(1,poissonNext)
        ciRatio <- ci(currentTime,lastTime)/poissonNext
        if (ciRatio > 1) stop("Problem with thinning.")
      }
    } else {
      currentTime <- currentTime + rexp(1,poissonRest)
      ciRatio <- ciRest(currentTime)/poissonRest
      if (ciRatio > 1) stop("Problem with thinning.")
      while (runif(1) > ciRatio) {
        currentTime <- currentTime + rexp(1,poissonRest)
        ciRatio <- ciRest(currentTime)/poissonRest
        if (ciRatio > 1) stop("Problem with thinning.")
      }
    }
    currentTime
  }

  while(currentTime <= duration) {
    currentTime <- nextTime(currentTime,lastTime)
    result[eventIdx] <- currentTime
    lastTime <- currentTime
    eventIdx <- eventIdx+1
    if (eventIdx > vLength) {
      result <- c(result,numeric(vLength))
      vLength <- length(result)
    }
  }
  result[result > 0]
  
}

## set the rng seed
set.seed(11006,"Mersenne-Twister")
## generate a log-normal pre train
preTrain <- cumsum(rlnorm(1000,-2.4,0.4))
preTrain <- preTrain[preTrain < 60]
## generate a post synaptic train with an
## exponentially decaying pre-synaptic excitation
post1 <- mkPostTrain(preTrain)
## generate a post synaptic train with a
## pulse-like pre-synaptic excitation
post2 <- mkPostTrain(preTrain,preFF=pulseFF)
## generate a post synaptic train with a
## pulse-like pre-synaptic inhibition
post3 <- mkPostTrain(preTrain,preFF=pulseFF,beta=-log(5))
## make a list of spikeTrain objects out of that
interData <- list(pre=as.spikeTrain(preTrain),
                  post1=as.spikeTrain(post1),
                  post2=as.spikeTrain(post2),
                  post3=as.spikeTrain(post3))
## remove the trains
rm(preTrain,post1,post2,post3)
## look at them
interData[["pre"]]
interData[["post1"]]
interData[["post2"]]
interData[["post3"]]
## compute cross-correlograms
interData.lt1 <- lockedTrain(interData[["pre"]],interData[["post1"]],laglim=c(-0.03,0.05),c(0,60))
interData.lt2 <- lockedTrain(interData[["pre"]],interData[["post2"]],laglim=c(-0.03,0.05),c(0,60))
interData.lt3 <- lockedTrain(interData[["pre"]],interData[["post3"]],laglim=c(-0.03,0.05),c(0,60))
## look at the cross-raster plots
interData.lt1
interData.lt2
interData.lt3
## look at the corresponding histograms
hist(interData.lt1,bw=0.0025)
hist(interData.lt2,bw=0.0025)
hist(interData.lt3,bw=0.0025)
## check out what goes on between post2 and post1
interData.lt1v2 <- lockedTrain(interData[["post2"]],interData[["post1"]],laglim=c(-0.03,0.05),c(0,60))
interData.lt1v2
hist(interData.lt1v2,bw=0.0025)

## fine
## create a GLM data frame using a 1 ms bin width
dfAll <- mkGLMdf(interData,delta=0.001,lwr=0,upr=60)
## build the sub-list relating to neuron 2
dfN2 <- dfAll[dfAll$neuron=="2",]
## fit dfN2 with a smooth effect for the elasped time since the last
## event of neuron 2 and another one with the elasped time since the
## last event from neuron 1. Use moroever only the events for which the
## the last event from neuron 1 occurred at most 100 ms ago.
dfN2.fit0 <- gam(event ~ s(lN.1,bs="cr") + s(lN.2,bs="cr"), data=dfN2, family=poisson, subset=(dfN2$lN.1 <=0.1))
## look at the summary
summary(dfN2.fit0)
## plot the smooth term of neuron 1
plot(dfN2.fit0,select=1,rug=FALSE,ylim=c(-0.8,0.8))
## Can you see the exponential presynatic effect with
## a 15 ms decay time appearing?
## Now check the dependence on lN.2
xx <- seq(0.001,0.3,0.001)
## plot the estimated conditional intensity when the last spike
## from neuron 1 came a long time ago (100 ms)
plot(xx,exp(predict(dfN2.fit0,data.frame(lN.1=rep(100,300)*0.001,lN.2=(1:300)*0.001))),type="l")
## add a line for the true conditional intensity
lines(xx,hlnorm(xx,-2.4,0.4)*0.001,col=2)
## do the same thing for the survival function
plot(xx,exp(-cumsum(exp(predict(dfN2.fit0,data.frame(lN.1=rep(100,300)*0.001,lN.2=(1:300)*0.001))))),type="l")
lines(xx,plnorm(xx,-2.4,0.4,lower.tail=FALSE),col=2)

## use gssanova
## split the data set in 2 parts, one for the fit, the other for the
## test
dfN2e <- dfN2[dfN2$time <= 20,]
dfN2l <- dfN2[dfN2$time > 20,]
## fit the same model as before with gssanova
dfN2.fit1 <- gssanova(event ~ lN.1 + lN.2, data=dfN2e, family="poisson", seed=20061001) 
## plot the effect of neuron 1
pred1 <- predict(dfN2.fit1,data.frame(lN.1=seq(0.001,0.220,0.001),
                                      lN.2=rep(median(dfN2e$lN.2),220)),
                 se=TRUE)
plot(seq(0.001,0.220,0.001),
     pred1$fit,type="l",
     ylim=c(min(pred1$fit-1.96*pred1$se.fit),max(pred1$fit+1.96*pred1$se.fit))
    )
lines(seq(0.001,0.220,0.001),pred1$fit-1.96*pred1$se.fit,lty=2)
lines(seq(0.001,0.220,0.001),pred1$fit+1.96*pred1$se.fit,lty=2)
## transform the time of the late part of the train
## first make sure than lN.1 and lN.2 are within the right bounds
m1 <- max(dfN2e$lN.1)
m2 <- max(dfN2e$lN.2)
dfN2l$lN.1 <- sapply(dfN2l$lN.1, function(x) min(m1,x))
dfN2l$lN.2 <- sapply(dfN2l$lN.2, function(x) min(m2,x))
predl <- predict(dfN2.fit1,dfN2l)
Lambda <- cumsum(exp(predl))
ttl <- mkCPSP(Lambda[dfN2l$event==1])
ttl
plot(summary(ttl))
## see what happens without time transformation
rtl <- mkCPSP(dfN2l$time[dfN2l$event==1])
plot(summary(rtl))

## Now repeat the fit including a possible contribution from neuron 3
dfN2.fit1 <- gam(event ~ s(lN.1,bs="cr") + s(lN.2,bs="cr") + s(lN.3,bs="cr"), data=dfN2, family=poisson, subset=(dfN2$lN.1 <=0.1) & (dfN2$lN.3 <= 0.1)) 
## Use the summary to see if the new element brings something
summary(dfN2.fit1)
## It does not!
## Now look at neurons 3 and 4 (ie, post2 and post3)
dfN3 <- dfAll[dfAll$neuron=="3",]
dfN3.fit0 <- gam(event ~ s(lN.1,k=20,bs="cr") + s(lN.3,k=15,bs="cr"),data=dfN3,family=poisson, subset=(dfN3$lN.1 <=0.1))
summary(dfN3.fit0)
plot(dfN3.fit0,select=1,ylim=c(-1.5,1.8),rug=FALSE)
dfN4 <- dfAll[dfAll$neuron=="4",]
dfN4.fit0 <- gam(event ~ s(lN.1,k=20,bs="cr") + s(lN.4,k=15,bs="cr"),data=dfN4,family=poisson, subset=(dfN4$lN.1 <=0.1))
summary(dfN4.fit0)
plot(dfN4.fit0,select=1,ylim=c(-1.8,1.5),rug=FALSE)
}
}
\keyword{ts}
