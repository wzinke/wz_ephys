\name{df4counts}
\alias{df4counts}
\title{Generates a Data Frame from a repeatedTrain Object After Time Binning}
\description{
  Generates a \code{\link{data.frame}} object out of a
  \code{repeatedTrain} object after time binning in order to study
  trials stationarity with a \code{\link{glm}} fit.
}
\usage{
df4counts(repeatedTrain, breaks = length(repeatedTrain))
}
\arguments{
  \item{repeatedTrain}{a \code{repeatedTrain} object or a list which can be
    coerced to such an object.}
  \item{breaks}{a numeric. A single number is interpreted has the number
    of bins; a vector is interpreted as the
    position of the "breaks" between bins.}
}
\details{
  The bins are placed between the \code{\link{floor}} of the smallest
  spike time and the \code{\link{ceiling}} of the largest one when
  \code{breaks} is a scalar. After time binning the number of spikes of
  each trial falling in each bin is counted (in the same way as the
  \code{counts} component of a \code{\link{psth}} list is
  obtained). This matrix of count is then formatted as a data frame.
}
\value{
  A \code{\link{data.frame}} with the following variables:
  \item{Count }{a count (number of spikes in a given bin at a given trial).}
  \item{Bin }{the bin index (a \code{\link{factor}}.}
  \item{Trial }{the trial index (a \code{\link{factor}}.}
  \item{Rate }{the count divided by the length of the corresponding
    bin.}
  \item{Time }{the time of the midpoints of the bins.}
}
\author{ Christophe Pouzat  \email{christophe.pouzat@gmail.com} }
\note{
  When a \code{\link{glm}} of the poisson family is used for subsequent
  analysis the important implicit hypothesis of an inhomogenous Poisson
  train is of course made.
}
\seealso{
  \code{\link{as.repeatedTrain}},
  \code{\link{psth}}
}
\examples{
\dontrun{
## Load the Vanillin responses of the first
## cockroach data set
data(CAL1V)
## convert them into repeatedTrain objects
## The stimulus command is on between 4.49 s and 4.99s
CAL1V <- lapply(CAL1V,as.repeatedTrain)
## Generate raster plot for neuron 1
raster(CAL1V[["neuron 1"]],c(4.49,4.99))
## make a smooth PSTH of these data
psth(CAL1V[["neuron 1"]],stimTimeCourse=c(4.49,4.99),breaks=c(bw=0.5,step=0.05),colCI=2,xlim=c(0,10))
## add a grid to the plot
grid()
## The response starts after 4.5 s and is mostly over after 6 s: create
## breaks accordingly
myBreaks <- c(0,2.25,4.5,seq(4.75,6.25,0.25),seq(6.5,11,0.5))
## get a count data frame
CAL1Vn1DF <- df4counts(CAL1V[["neuron 1"]],myBreaks)
## use a box plot to look at the result
boxplot(Rate ~ Time, data=CAL1Vn1DF)
## watch out here the time scale is distorted because of our
## choice of unequal bins
## Fit a glm of the Poisson family taking both Bin and Trial effects
CAL1Vn1DFglm <- glm(Count ~ Bin + Trial,family=poisson,data=CAL1Vn1DF)
## use an anova to see that both the Bin effect and the trial effect are
## highly significant
anova(CAL1Vn1DFglm, test="Chisq")
}}
\keyword{ts}
