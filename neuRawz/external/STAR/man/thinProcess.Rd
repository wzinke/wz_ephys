\name{thinProcess}
\Rdversion{1.1}
\alias{thinProcess}
\alias{maxIntensity}
\alias{mkSelf}
\alias{mkMappedI}
\alias{mkSimFct}
\alias{mkPostSimAnalysis}

\title{
Simulate and Analyse Data From a Model Fitted With gss
}
\description{
Functions \code{thinProcess} simulates a spike train and \code{mkSimFct}
returns a simulating function from a \code{gssanova} fitted
model. Ogata's thinning simulation 
method is used. Functions \code{maxIntensity},
\code{mkSelf} and \code{mkMappedI} are utility functions for the first
two. Function \code{mkPostSimAnalysis} returns a function analysing a
simulated spike train. Functions \code{mkSimFct} and
\code{mkPostSimAnalysis} return functions which can in principle be
safely used in parallel applications, that is, they have everything they
need in their closure. 
}
\usage{
thinProcess(object, m2uFctList, trueData, formerSpikes,
            intensityMax, ...)
maxIntensity(object, dfWithTime, ...)
mkSelf(m2uSelf)
mkMappedI(m2uI, lag = 1)
mkSimFct(object, m2uFctList, trueData, formerSpikes,
         intensityMax, ...)
mkPostSimAnalysis(stList, trainNumber = 1, timeWindow,
                  objects, dfFct)
}
\arguments{
  \item{object}{A \code{\link[gss]{ssanova}} or a
    \code{\link[gss]{ssanova0}} object.} 
  \item{m2uFctList}{A list of functions. There should be as many functions as
    there are "internal" variables in \code{object}. An internal
    variable is a variable 
    whose value is changed by the occurrence of a spike, like the elapsed
    time since the last spike of the simulated neuron, the duration of a
    former inter spike interval of a given lag, etc. The names of the
    components (functions) of the list should be the names of the variables.
    Each function should correspond to the map to uniform function used before
    fitting the model.}
  \item{m2uSelf}{The map to uniform function used to transform the
    actual elapsed time since the last spike values before fitting the model.}
  \item{m2uI}{The map to uniform function used to transform the actual
    former ISI durations before fitting the model.}
  \item{lag}{The considered lag (integer > 0).}
  \item{trueData}{A data frame containing the "true data" of the
    simulated epoch. This is to ensure that "external" variables such as
    the elapsed time since the last spike of a functionally coupled
    neuron are available.}
  \item{formerSpikes}{A vector of previous spike times. This is to make the
    computation of former inter spike intervals possible in every case.}
  \item{intensityMax}{The value of the maximal intensity. If missing function
    \code{maxIntensity} is called to estimate it.}
  \item{dfWithTime}{A data frame with one variable named "time". The
    latter variable is used to obtain the bin width with which the
    original spike train was discretized.}
  \item{stList}{The list of \code{spikeTrain} objects with one of the
    trains partly simulated. A single (partly simulated)
    \code{spikeTrain} object can also be used.}
  \item{trainNumber}{An integer, the index of the modeled and simulated
    spike train in \code{stList}.}
  \item{timeWindow}{A numeric vector of length 1 or 2. This argument
    specifies the time domain over which the fits contained in argument
    objects was performed. It is implicitly assumed that the (partial)
    simulation was performed outside this time domain. When a vector of
    length 1 is used the fitting time domain is taken as
    \code{c(0,timeWindow)}.}
  \item{objects}{A list of \code{\link[gss]{ssanova}} or
    \code{\link[gss]{ssanova0}} objects. Each element of
    the list is a "model" with which analysis will be performed. A
    single \code{\link[gss]{ssanova}} or \code{\link[gss]{ssanova0}}
    object can also be used.}
  \item{dfFct}{A function whose argument is a the same as the first
    argument of function \code{\link[STAR]{mkGLMdf}} and which returns a
    data frame suitable for use as argument \code{newdata} in
    \code{\link[gss]{predict.ssanova}}.}
  \item{\dots}{Additional arguments passed to \code{\link{optim}} which
    is called internally with the "BFGS" method in function
    \code{maxIntensity}. In functions \code{thinProcess} and
    \code{mkSimFct}, additional argument passed to function
    \code{maxIntensity} if necessary.}
}
\details{
  Function \code{thinProcess} simulates a spike train with Ogata's
  thinning method (Ogata, 1981). The latter method required the maximal
  intensity of the process to be known. If such is not the case, that
  is, if argument \code{intensityMax} is left \code{missing}, a proposed
  maximal intensity is obtained with function \code{maxIntensity}. If
  during the simulation an actual intensity larger than the given
  \code{intensityMax} occurs, the simulation is interrupted and an error
  message is generated.

  Function \code{maxIntensity} uses the central point of the variable
  space as its intial guess. The "BFGS" method of \code{\link{optim}} is
  used to find the maximal intensity.
  
  Function \code{mkPostSimAnalysis} uses function
  \code{\link[codetools]{findGlobals}} in order to find among the
  functions called by \code{dfFct} the ones which are defined in the
  global environment. These functions are copied in the environment
  (Gentleman and Ihaka, 2000) of the function returned by
  \code{mkPostSimAnalysis}. If the global environment defined function
  called by \code{dfFct} \emph{do not} call themselves over functions
  defined in the global environment, the returned function can be safely
  used as argument \code{fun} of package \code{snow}'s \code{clusterApply}
  function.    
}
\value{
  \code{thinProcess} returns a \code{spikeTrain} object.

  \code{maxIntensity} returns the "proposed" maximal intensity (in Hz).

  \code{mkSelf} returns a \code{function} taking two arguments:
  \code{self(proposedtime,st)}.

  \code{mkMappedI} returns a \code{function} taking two arguments:
  \code{function(proposedtime,st)}.

  \code{mkSimFct} returns a function simulating a \code{spikeTrain}
  object. The simulation is done with function \code{thinProcess}. The
  returned function takes no argument. The maximal intensity required by
  the thinning method is stored in the closure of the returned function.

  \code{mkPostSimAnalysis} returns a function taking a \code{spikeTrain}
  object as its single argument. This function returns a list of
  lists. Each list correspond to one of the models in argument
  \code{objects}. Each sub list has two components: \code{lpp} (the log
  predictive probability) and \code{ttt} (the time transformed train, a
  CountingProcessSamplePath object).
}
\references{
  Gentleman, R. and Ihaka, R. (2000) Lexical Scope and Statistical
  Computing. \emph{Journal of Computational and Graphical Statistics}
  \bold{9}: 491-508.
  
  Ogata, Y. (1981) On Lewis' simulation method for point
  processes. \emph{IEEE Transactions on Information Theory}
  \bold{IT-29}: 23-31. 
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com}}

\note{These functions are designed to implement a rather specific type
  of analysis which is exposed in the "big STAR tutorial" available at:
  \url{http://sites.google.com/site/spiketrainanalysiswithr/}. The
  exemple below shows a "complete" analysis, more details and other
  exemples can be found in the big tutorial.
}
\seealso{
  \code{\link[gss]{gssanova}},
  \code{\link{as.spikeTrain}},
  \code{\link{mkGLMdf}},
  \code{\link{mkCPSP}},
}
\examples{
\dontrun{
## The run times given in the sequel were measured on a laptop
## with a dual core CPU: 2x Intel Core 2 Duo CPU P9500 @ 2.53GHz
## The RAM was 4 GB large. The PC ran Ubuntu 9.04 and R-2.9.2
## compiled with a dynamically link ATLAS version of BLAS.

## Start by loading the data set into the work space.
data(e060824spont)
## Get a summary of neuron 1 spike train. 
summary(e060824spont[["neuron 1"]])
## Run an automatic analysis of the train (takes ~ 4.22 s)
reportHTML(e060824spont[["neuron 1"]],filename="e060824spont_1",otherST=e060824spont[c(2)],maxiter=100)
## The renewal tests show that the discharge is not the one of
## a renewal process. The cross-correlogram shows no sign of
## coupling between the two neurons of the data set.

## Compute the partial autocorrelation function to get an idea
## of how many previous interspike intervals (ISIs) should be
## included in the model.
acf(diff(e060824spont[["neuron 1"]]),type="partial")
## The pacf plot suggests that the last ISI should be enough.
## Build the data frame.
DFA <- mkGLMdf(e060824spont[["neuron 1"]],0.004,0,59)
DFA <- within(DFA,i1 <- isi(DFA,lag=1))
DFA <- DFA[complete.cases(DFA),]

## look a the ECDF of the elapsed time since the last spike,
## that is, variable lN.1 of the data frame and of the last
## ISI (variable i1 of the data frame).
layout(matrix(1:2,nc=2))
with(DFA,plot(ecdf(lN.1),pch="."))
with(DFA,lines(range(lN.1),c(0,1),col=2,lty=2))
with(DFA,plot(ecdf(i1),pch="."))
with(DFA,lines(range(i1),c(0,1),col=2,lty=2))
## The distributions of these varaibles are clearly (and not
## surprisingly) non-uniform.

## Build emprirical functions mapping the two variables to uniform
## ones
m2u1 <- mkM2U(DFA,"lN.1",0,28.5)
m2ui <- mkM2U(DFA,"i1",0,28.5,maxiter=200)
DFA <- within(DFA,e1t <- m2u1(lN.1))
DFA <- within(DFA,i1t <- m2ui(i1))
## Cehck that the transformations did their job.
with(DFA,plot(ecdf(e1t),pch="."))
with(DFA,lines(range(e1t),c(0,1),col=2,lty=2))
with(DFA,plot(ecdf(i1t),pch="."))
with(DFA,lines(range(i1t),c(0,1),col=2,lty=2))

## The heavy computations to follow will be performed
## in parallel using the snow package.
library(snow)
## Define the number of slaves
nbSlaves <- 2
## Create the cluster.
cl <- makeCluster(rep("localhost",nbSlaves),type="SOCK")
## load STAR on each slave.
clusterEvalQ(cl,library(STAR))

## Define a function making a function performing the 
## fit with gssanova and suitable for a parallel implementation.
## The returned function does in addition time transform the
## spike train on the time window not used for the fit.
mkPFct <- function(df=DFA) {
  df
  PFct <- function(gtime,
                   fmla=event~e1t*i1t,
                   seed=20061001) {
    GF <- gssanova(fmla,
                   data=subset(df, time \%in\% gtime),
                   family="binomial",
                   seed=seed)
    tt <- GF \%tt\% subset(df, !(time \%in\% gtime))
    list(GF=GF,tt=tt)
  }
  PFct
}

PFct1 <- mkPFct()
## Now PFct1 returns a list with two elements: the "fit object" (GF)
## and the time transformed train (tt)

## Create a list suitable as the second argument for clusterApply
gtList <- list(early=with(DFA,time[time<=29.5]),
               late=with(DFA,time[time>29.5])
               )

## Fit and test a model with interaction between the (mapped)
## ellasped time since the last spike and the (mapped) last
## ISI. This takes ~ 95 s.
GF1.e060824spont.1 <- clusterApply(cl, gtList, PFct1)

## Look a the test battery
plot(summary(GF1.e060824spont.1[[1]][[2]]),which=c(1,2,4,6))
plot(summary(GF1.e060824spont.1[[2]][[2]]),which=c(1,2,4,6))

## Fit and test a model without interaction between the (mapped)
## ellasped time since the last spike and the (mapped) last
## ISI. This takes ~ 61 s.
GF2.e060824spont.1 <- clusterApply(cl, gtList, PFct1,fmla=event ~ e1t+i1t)

## Look a the test battery
plot(summary(GF2.e060824spont.1[[1]][[2]]),which=c(1,2,4,6))
plot(summary(GF2.e060824spont.1[[2]][[2]]),which=c(1,2,4,6))

## Compute the "predictive log probability" with Model 2
## (without interaction). This takes ~ 1.6 s
(GF2.e060824spont.1.logProb <- predictLogProb(GF2.e060824spont.1[[1]][[1]],subset(DFA,time>29.5))+predictLogProb(GF2.e060824spont.1[[2]][[1]],subset(DFA,time<=29.5)))

## Compute the "predictive log probability" with Model 1
## (with interaction). This takes ~ 3.5 s
(GF1.e060824spont.1.logProb <- predictLogProb(GF1.e060824spont.1[[1]][[1]],subset(DFA,time>29.5))+predictLogProb(GF1.e060824spont.1[[2]][[1]],subset(DFA,time<=29.5)))

## Prepare the simulations using Model 1 and 2
## Define a function initializing a mrg32k3a RNG from
## the rstream package on each slave
initMRG32k3a <- function(cl) {
  clusterEvalQ(cl,library(rstream))
  invisible(clusterCall(cl,
                        function() {
                          cmd <- parse(text=".s <- new(\"rstream.mrg32k3a\")")
                          eval(cmd,env=globalenv())
                        }
                        )
            )
  cat(paste(paste(clusterEvalQ(cl,rstream.sample(.s)),collapse=","),"\n"))
  invisible(clusterEvalQ(cl,rstream.reset(.s)))
}

## Define a function returning a list of independent and packed
## mrg32k3a rngs.
mkLecuyerList <- function(cl, ## a snow cluster
                          seed,
                          ...) {
  nbWorkers <- length(cl)
  lecuyerList <- vector(mode="list",length=nbWorkers)
  for (wIdx in 1:nbWorkers) {
    if (wIdx == 1) {
      if (!missing(seed)) lecuyerList[[1]] <- new("rstream.mrg32k3a",seed=seed)
      else lecuyerList[[1]] <- new("rstream.mrg32k3a")
    } else lecuyerList[[wIdx]] <- new("rstream.mrg32k3a")
    rstream.packed(lecuyerList[[wIdx]]) <- TRUE
  }
  lecuyerList
}

## Define a function setting the uniform rng of each slave
## to one of the independent mrg32k3a rngs created by
## mkLecuyerList.
clusterSetupRSTREAM <- function(cl, ## a snow cluster
                                lecuyerList
                                ) {


  setLecuyer <- function(packedlecuyer) {
    assign("lecuyer",packedlecuyer,env=globalenv())
    cmd <- parse(text="rstream.packed(lecuyer)<-FALSE")
    eval(cmd,env=globalenv())
  }

  clusterApply(cl,lecuyerList,setLecuyer)
  clusterEvalQ(cl,rstream.RNG(lecuyer))

}

## Load package rstream on master.
library(rstream)
## Initialize mrg32k3a rngs on each slave.
initMRG32k3a(cl)
## Create the list of independent mrg32k3a rngs
## on master.
theList <- mkLecuyerList(cl,seed=rep(20061001,6))
## Set the uniform rng of each slave to one of the
## independent mrg32k3a rngs just created.
clusterSetupRSTREAM(cl,theList)

## Define a list of map to uniform functions
fList.e060824spont.1 <- list(e1t=mkSelf(m2u1), i1t=mkMappedI(m2ui))
## Define a simulating function from Model 1 fitted on the
## half of the data set.
simF1.e060824spont.1 <- mkSimFct(object=GF1.e060824spont.1[[1]][[1]],
                                 m2uFctList=fList.e060824spont.1,
                                 trueData=subset(DFA,time>29.5),
                                 formerSpikes=with(DFA,time[event==1][time[event==1] <= 29.5])
                                 )

## Define a simulating function from Model 2 fitted on the
## half of the data set.
simF2.e060824spont.1 <- mkSimFct(object=GF2.e060824spont.1[[1]][[1]],
                                 m2uFctList=fList.e060824spont.1,
                                 trueData=subset(DFA,time>29.5),
                                 formerSpikes=with(DFA,time[event==1][time[event==1] <= 29.5])
                                 )

## Define the number of simulations to carry out.
nbRep <- 100
## Simulate spike trains in parallel using Model 1.
## This takes ~ 670 s.
sim1.e060824spont.1 <- clusterApply(cl,
                                    rep(nbRep/nbSlaves,nbSlaves),
                                    function(n,SF) lapply(1:n, function(idx) SF()),
                                    SF=simF1.e060824spont.1)

## Convert the returned list of lists into a single
## big list.
sim1.e060824spont.1 <- c(sim1.e060824spont.1[[1]],
                         sim1.e060824spont.1[[2]])

## Simulate spike trains in parallel using Model 1.
## This takes ~ 425 s.
sim2.e060824spont.1 <- clusterApply(cl,
                                    rep(nbRep/nbSlaves,nbSlaves),
                                    function(n,SF) lapply(1:n, function(idx) SF()),
                                    SF=simF2.e060824spont.1)

## Convert the returned list of lists into a single
## big list.
sim2.e060824spont.1 <- c(sim2.e060824spont.1[[1]],
                         sim2.e060824spont.1[[2]])

## Define a function generating automatically the
## proper data frame from the simulated data.
mkDF.e060824spont.1 <- function(stList) {
  DF <- mkGLMdf(stList,0.004,0,59)
  DF <- within(DF,i1 <- isi(DF,lag=1))
  DF <- DF[complete.cases(DF),]
  DF <- within(DF,e1t <- m2u1(lN.1))
  DF <- within(DF,i1t <- m2ui(i1))
  DF
}

## Define a function analysis the simulated trains with
## both Model 1 and 2.
PSAFct <- mkPostSimAnalysis(e060824spont[[1]],1,29.5,list(GF1.e060824spont.1[[1]][[1]],GF2.e060824spont.1[[1]][[1]]),mkDF.e060824spont.1)

## Analyze the simulations done with Model 1.
## This takes ~ 400 s
sim1.e060824spont.1.psa <- clusterApply(cl,sim1.e060824spont.1,PSAFct)

## Analyze the simulations done with Model 2.
## This takes ~ 400 s
sim2.e060824spont.1.psa <- clusterApply(cl,sim2.e060824spont.1,PSAFct)

## Get the log predictive probability assuming Model 1 for
## simulations done with Model 1.
sim1.e060824spont.1.lpp1 <- sapply(sim1.e060824spont.1.psa, function(l) l[[1]]$lpp)
## Get the log predictive probability assuming Model 2 for
## simulations done with Model 1.
sim1.e060824spont.1.lpp2 <- sapply(sim1.e060824spont.1.psa, function(l) l[[2]]$lpp)
## Get the log predictive probability assuming Model 1 for
## simulations done with Model 2.
sim2.e060824spont.1.lpp1 <- sapply(sim2.e060824spont.1.psa, function(l) l[[1]]$lpp)
## Get the log predictive probability assuming Model 2 for
## simulations done with Model 2.
sim2.e060824spont.1.lpp2 <- sapply(sim2.e060824spont.1.psa, function(l) l[[2]]$lpp)

## Get the observed log predictive probability with each model.
e060824spont.1.lpp1 <- predictLogProb(GF1.e060824spont.1[[1]][[1]],subset(DFA,time>29.5))
e060824spont.1.lpp2 <- predictLogProb(GF2.e060824spont.1[[1]][[1]],subset(DFA,time>29.5))

## Get the difference of observed log predictive probabilities.
e060824spont.1.lppDiff <- e060824spont.1.lpp1 - e060824spont.1.lpp2

## Look at the correlation between the log predictive probabilities
## obtained with Model 1 and 2 with data simulated with Model 1.
plot(sim1.e060824spont.1.lpp1,sim1.e060824spont.1.lpp2,main="log prob with M2 vs log prob with M1 when M1 is true",xlab="log prob with M1",ylab="log prob with M2")

## Plot the ECDF of the log predictive probabilities obtained
## with Model 1 with data simulated with Model 1.
plot(ecdf(sim1.e060824spont.1.lpp1),pch=".",main="log prob with Model 1 when Model 1 is true")
## Show the observed value of this statistic.
segments(e060824spont.1.lpp1,0,e060824spont.1.lpp1,sum(sim1.e060824spont.1.lpp1 <= e060824spont.1.lpp1)/nbRep,col=2,lwd=2)
segments(-1600,sum(sim1.e060824spont.1.lpp1 <= e060824spont.1.lpp1)/nbRep,e060824spont.1.lpp1,sum(sim1.e060824spont.1.lpp1 <= e060824spont.1.lpp1)/nbRep,col=2,lwd=2)
## Plot the ECDF of the log predictive probabilities obtained
## with Model 2 with data simulated with Model 1.
plot(ecdf(sim1.e060824spont.1.lpp2),pch=".",main="log prob with Model 2 when Model 1 is true")
## Show the observed value of this statistic.
segments(e060824spont.1.lpp2,0,e060824spont.1.lpp2,sum(sim1.e060824spont.1.lpp2 <= e060824spont.1.lpp2)/nbRep,col=2,lwd=2)
segments(-1800,sum(sim1.e060824spont.1.lpp2 <= e060824spont.1.lpp2)/nbRep,e060824spont.1.lpp2,sum(sim1.e060824spont.1.lpp2 <= e060824spont.1.lpp2)/nbRep,col=2,lwd=2)

## Plot the ECDF of the difference of the log predictive probabilities 
## obtained with data simulated with Model 1.
plot(ecdf(sim1.e060824spont.1.lpp1-sim1.e060824spont.1.lpp2),pch=".",main="log prob with M1 - log prob with M2 when M1 is true")
## Show the observed value of this statistic.
segments(e060824spont.1.lppDiff,0,e060824spont.1.lppDiff,sum(sim1.e060824spont.1.lpp1-sim1.e060824spont.1.lpp2<=e060824spont.1.lppDiff)/nbRep,col=2,lwd=2)
segments(-10,sum(sim1.e060824spont.1.lpp1-sim1.e060824spont.1.lpp2<=e060824spont.1.lppDiff)/nbRep,e060824spont.1.lppDiff,sum(sim1.e060824spont.1.lpp1-sim1.e060824spont.1.lpp2<=e060824spont.1.lppDiff)/nbRep,col=2,lwd=2)

## Look at the correlation between the log predictive probabilities
## obtained with Model 1 and 2 with data simulated with Model 2.
plot(sim2.e060824spont.1.lpp1,sim2.e060824spont.1.lpp2,main="log prob with M2 vs log prob with M1 when M2 is true",xlab="log prob with M1",ylab="log prob with M2")
## Plot the ECDF of the log predictive probabilities obtained
## with Model 1 with data simulated with Model 2.
plot(ecdf(sim2.e060824spont.1.lpp1),pch=".",main="log prob with Model 1 when Model 2 is true")
## Show the observed value of this statistic.
segments(e060824spont.1.lpp1,0,e060824spont.1.lpp1,sum(sim2.e060824spont.1.lpp1 <= e060824spont.1.lpp1)/nbRep,col=2,lwd=2)
segments(-2000,sum(sim2.e060824spont.1.lpp1 <= e060824spont.1.lpp1)/nbRep,e060824spont.1.lpp1,sum(sim2.e060824spont.1.lpp1 <= e060824spont.1.lpp1)/nbRep,col=2,lwd=2)
## Plot the ECDF of the log predictive probabilities obtained
## with Model 2 with data simulated with Model 2.
plot(ecdf(sim2.e060824spont.1.lpp2),pch=".",main="log prob with Model 2 when Model 2 is true")
## Show the observed value of this statistic.
segments(e060824spont.1.lpp2,0,e060824spont.1.lpp2,sum(sim2.e060824spont.1.lpp2 <= e060824spont.1.lpp2)/nbRep,col=2,lwd=2)
segments(-2000,sum(sim2.e060824spont.1.lpp2 <= e060824spont.1.lpp2)/nbRep,e060824spont.1.lpp2,sum(sim2.e060824spont.1.lpp2 <= e060824spont.1.lpp2)/nbRep,col=2,lwd=2)
## Plot the ECDF of the difference of the log predictive probabilities 
## obtained with data simulated with Model 1.
## Make sure that the scale of the x axis is right.
xlim=c(min(c(-e060824spont.1.lppDiff,sim2.e060824spont.1.lpp2-sim2.e060824spont.1.lpp1)),max(sim2.e060824spont.1.lpp2-sim2.e060824spont.1.lpp1))
plot(ecdf(sim2.e060824spont.1.lpp2-sim2.e060824spont.1.lpp1),pch=".",main="log prob with M2 - log prob with M1 when M2 is true",xlim=xlim)
## Show the observed value of this statistic.
points(-e060824spont.1.lppDiff,0,pch=16,col=2)

## Stop the snow cluster.
stopCluster(cl)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
\keyword{htest}
\keyword{datagen}
