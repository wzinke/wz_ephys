\name{frt}
\alias{frt}
\alias{\%frt\%}
\title{Computes Forward Recurrence Times from Two transformedTrain Objects}
\description{
  Computes the (transformed) time differences between spikes of a
  \code{refTrain} and the (next) ones of a \code{testTrain}. Both
  \code{refTrain} and \code{testTrain} should be \code{transformedTrain} objects.
  
}
\usage{
frt(refTrain, testTrain)
refTrain \%frt\% testTrain
}
\arguments{
  \item{refTrain}{a \code{transformedTrain} object.}
  \item{testTrain}{a \code{transformedTrain} object.}
}
\details{
  When two spike trains have been time transformed using \emph{the same}
  procedure, which does make one of the trains (the \code{testTrain}) the realization a
  homogeneous Poisson process with rate 1, the elapsed time between the
  spikes of the other train (\code{refTrain}) and the ones of
  \code{testTrain} should be exponentially distributed with rate
  1. These elapsed times are returned by \code{frt}.
}
\value{
  An object of class \code{frt} containing the elapsed times.
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com}}
\seealso{
  \code{\link{transformedTrain}},
  \code{\link{plot.frt}},
  \code{\link{summary.frt}},
  \code{\link{mkGLMdf}}
}
\examples{
\dontrun{
## Let us consider neuron 1 of the CAL2S data set
data(CAL2S)
CAL2S <- lapply(CAL2S,as.spikeTrain)
CAL2S[["neuron 1"]]
renewalTestPlot(CAL2S[["neuron 1"]])
summary(CAL2S[["neuron 1"]])
## Make a data frame with a 4 ms time resolution
cal2Sdf <- mkGLMdf(CAL2S,0.004,0,60)
## keep the part relative to neuron 1, 2 and 3 separately
n1.cal2sDF <- cal2Sdf[cal2Sdf$neuron=="1",]
n2.cal2sDF <- cal2Sdf[cal2Sdf$neuron=="2",]
n3.cal2sDF <- cal2Sdf[cal2Sdf$neuron=="3",]
## remove unnecessary data
rm(cal2Sdf)
## Extract the elapsed time since the second to last and
## third to last for neuron 1. Normalise the result. 
n1.cal2sDF[c("rlN.1","rsN.1","rtN.1")] <- brt4df(n1.cal2sDF,"lN.1",2,c("rlN.1","rsN.1","rtN.1"))
## load mgcv library
library(mgcv)
## fit a model with a tensorial product involving the last
## three spikes and using a cubic spline basis for the last two
## To gain time use a fixed df regression spline
n1S.fitA <- gam(event ~ te(rlN.1,rsN.1,bs="cr",fx=TRUE) + rtN.1,data=n1.cal2sDF,family=binomial(link="logit"))
## transform time
N1.Lambda <- transformedTrain(n1S.fitA)
## check out the resulting spike train using the fact
## that transformedTrain objects inherit from spikeTrain
## objects
N1.Lambda
## Use more formal checks
summary(N1.Lambda)
plot(N1.Lambda,which=c(1,2,4,5),ask=FALSE)
## Transform spike trains of neuron 2 and 3
N2.Lambda <- transformedTrain(n1S.fitA,n2.cal2sDF$event)
N3.Lambda <- transformedTrain(n1S.fitA,n3.cal2sDF$event)
## Check interactions
summary(N2.Lambda \%frt\% N1.Lambda)
summary(N3.Lambda \%frt\% N1.Lambda)
plot(N2.Lambda \%frt\% N1.Lambda,ask=FALSE)
plot(N3.Lambda \%frt\% N1.Lambda,ask=FALSE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
