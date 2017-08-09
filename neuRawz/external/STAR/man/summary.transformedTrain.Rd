\name{summary.transformedTrain}
\alias{summary.transformedTrain}
\title{Summary of transformedTrain Objects}
\description{
  Generates a concise summary of \code{transformedTrain} objects. It is mostly
  intended for use in batch processing situations where a decision to stop with the
  current model or go on with a more complicated one must be made automatically.
}
\usage{
\method{summary}{transformedTrain}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a \code{\link{transformedTrain}} object.}
  \item{\dots}{additional arguments passed to \code{\link{varianceTime}}.}
}
\details{
  \code{summary.transformedTrain} computes summary statistics
  corresponding to plot 1, 2 and 5 of
  \code{\link{plot.transformedTrain}}.

  The first plot tests the uniformity of the spikes (transformed) times
  on the (transformed) observation window using a KS test. If the
  ecdf of the (transformed) times is within 
  the 95\% band then the first element of component \code{uniformOnTTime}
  of the returned list is set to \code{TRUE}. It is set to \code{FALSE}
  otherwise. The second component is relative to the 99\% band.

  The second plot tests the exponential distribution of the intervals
  between successive spikes transformed times. Again if the empirical
  curve stays within the 95, respectively 99\%, confidence band, the
  first, respectively second, element of component \code{BermanTest}
  of the returned list is set to \code{TRUE}. It is set to \code{FALSE}
  otherwise.

  The fifth plot tests that the variance is equal to the length of the
  (transformed) observation time for \code{object}, using point-wise
  CI. If n different observation times are defined over the whole
  observation window, we expect (1 - CI/100)*n points
  to be out with an approximate binomial distribution. For each CI
  defined (95 and 99\%, by default), component \code{VarTime}
  of the returned list contains the probability of observing a number as
  large as or smaller than the one observed under the binomial null hypothesis. 
  
}
\value{
  A \code{list} with the following 3 components:
  \item{uniformOnTTime }{A two named components vector of boolean.}
  \item{BermanTest }{A two named components vector of boolean.}
  \item{VarTime }{A named component vector with as many components as
    passed to \code{\link{varianceTime}} via the \code{...} argument
    with p-values of a binomial distribution.}
  
}
\references{
  Ogata, Yosihiko (1988) Statistical Models for Earthquake Occurrences and Residual
  Analysis for Point Processes. \emph{Journal of the American
    Statistical Association} \bold{83}: 9-27.

  Brown, E. N., Barbieri, R., Ventura, V., Kass, R. E. and Frank,
  L. M. (2002) The time-rescaling theorem and its application to neural
  spike train data analysis. \emph{Neural Computation} \bold{14}:
  325-346.
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com}}
\seealso{
  \code{\link{transformedTrain}},
  \code{\link{plot.transformedTrain}},
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
\keyword{models}
\keyword{smooth}
\keyword{regression}
