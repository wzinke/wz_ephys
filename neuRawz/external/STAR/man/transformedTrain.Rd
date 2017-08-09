\name{transformedTrain}
\alias{transformedTrain}
\title{Performs Time Transformation of Spike Trains Fitted with glm or gam}
\description{
  Transform spike times from a \code{glm} or \code{gam} fitted model as
  defined by Ogata (1988) and Brown et al (2002). If the model structure
  is "correct" and if the model parameters are properly estimated the
  result of the time transformation should be the realization of a
  Poisson process with rate 1.
}
\usage{
transformedTrain(obj, target = obj$data$event, select)
}
\arguments{
  \item{obj}{An object returned by \code{\link{gam}} or \code{\link{glm}}.}
  \item{target}{A binary (0,1) vector of integers with the same length
    as \code{dim(obj$data)[1]} or a vector of indexes giving the
    discretized times of events. All these indexes should then be
    included in \code{seq(dim(obj$data)[1])}.}
  \item{select}{A character string defining a condition to be fulfilled
  by the event in order to be selected, like: \code{time <= 6}. This is
  evaluated after parsing in the data frame of \code{obj}.}
}
\details{
  The \code{fitted.values} component of \code{obj} contains the
  (estimated) probability to observe a spike in each time bin where the covariates
  required by the fitted model were defined. It is then straightforward
  to show using the concept of \emph{product integral} (Kalbfleisch and
  Prentice, 2002; Andersen et al, 1993),provided that the time bin width
  is small enough to have a very small probability in each bin, that the
  cumulated sum of these probabilities is the expected number of events
  observed up to a given time. This expected number of events
  which is returned by \code{transformedTrain}. It is also the result of
  the "time transformation" proposed by Ogata (1988) and brought to the
  spike train analysis field under the name "time rescaling (theorem)"
  by Brown et al (2002).

  \code{transformedTrain} can also be used to transform the times of the
  spikes of neurons whose spike trains were simultaneously recorded and
  discretized \emph{in exactly the same way} as the neuron used to
  generate \code{obj}. This is useful to explore the possibility of
  functional interactions between a putative pre-synaptic neuron (whose
  spike train would correspond to argument \code{target}) and a
  post-synaptic one used to generate \code{obj}.
}
\value{
  \code{transformedTrain} returns an object of class
  \code{transformedTrain} inheriting from class \code{spikeTrain}. The
  object is fundamentally a numeric vector with strictly increasing
  elements containing the transformed times (or the expected number of events).
}
\references{
  Ogata, Yosihiko (1988) Statistical Models for Earthquake Occurrences and Residual
  Analysis for Point Processes. \emph{Journal of the American
    Statistical Association} \bold{83}: 9-27.

  Brown, E. N., Barbieri, R., Ventura, V., Kass, R. E. and Frank,
  L. M. (2002) The time-rescaling theorem and its application to neural
  spike train data analysis. \emph{Neural Computation} \bold{14}:
  325-346.

  Kalbfleisch, John D. and Prentice, Ross L. (2002) \emph{The
    Statistical Analysis of Failure Time Data.} Wiley Interscience.
  
  Andersen, Per Kragh, Borgan, Ornulf, Gill, Richard D. and Keiding,
  Niels (1993) \emph{Statistical Models Based on Counting Processes.}
  Springer-Verlag. 
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com}}
\note{
  As mentioned only the spikes for which the covariates of the model are
  available have their times transformed. That practically means that
  the length of the \code{transformedTrain} object returned by function
  \code{transformedTrain} \emph{can be shorter} than the length of the
  original \code{spikeTrain} object (or more precisely than the number
  of spikes defined in \code{target}). If one works with a model
  involving the elapsed times since the last three spikes then the
  fourth spike of the train will be the first to be transformed. You
  should therefore expect some left truncation of the data at the
  beginning of each acquisition epoch.
}
\seealso{
  \code{\link{plot.transformedTrain}},
  \code{\link{summary.transformedTrain}},
  \code{\link{mkGLMdf}},
  \code{\link{data.frame}},
  \code{\link{glm}},
  \code{\link{mgcv}}
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
