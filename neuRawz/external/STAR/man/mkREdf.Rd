\name{mkREdf}
\alias{mkREdf}
\title{Evaluates RateEvolutions for spikeTrain Lists and Returns Data Frame}
\description{
  Given a list of \code{spikeTrain} or \code{repeatedTrain} objects
  \code{mkREdf} evaluates the rate evolution of each train and returns a
  data frame suitable for use with \code{coplot}, \code{xyplot} and \code{qplot}.
}
\usage{
mkREdf(x, longitudinal, across, bw,
       kernel=c("gaussian", "epanechnikov", "rectangular",
                "triangular", "biweight", "cosine", "optcosine"),
       n=512, from, to, na.rm=FALSE, minusMean=FALSE)
}
\arguments{
  \item{x}{a \emph{named} list of \code{spikeTrain} or
    \code{repeatedTrain} objects.}
  \item{longitudinal}{a \code{character} vector with the names of the
    different "conditions" applied to each neuron like "ctl", "bicu" or
    "stim. 1", "stim. 2", ..., "stim. 20". Default provided.}
  \item{across}{a \code{character} vector with the names of the
    different neurons. Default provided.}
  \item{bw}{see \code{\link{rateEvolution}}. This can be a vector.}
  \item{kernel}{see \code{\link{rateEvolution}}.}
  \item{n}{see \code{\link{rateEvolution}}.}
  \item{from}{see \code{\link{rateEvolution}}.}
  \item{to}{see \code{\link{rateEvolution}}.}
  \item{na.rm}{see \code{\link{rateEvolution}}.}
  \item{minusMean}{ should the mean of the rate evolution along the
    across "dimension" be subtracted from each individual rate evolution
    along this dimension?}
}
\details{
  \code{mkREdf} calls \code{\link{rateEvolution}} on every
  \code{spikeTrain} in \code{x}. If \code{from} and \code{to} are
  missing, they are internally set to the \code{floor} of the global
  minimal spike time contained in \code{x} and to the \code{ceiling} of
  the global maximal time.
}
\value{
  A data frame with the following variables:
  \item{time }{The time (in s) at which the rate was evaluated.}
  \item{rate }{The rate (in 1/s).}
  \item{longitudinal }{A factor corresponding to the argument with the
    same name.}
  \item{across }{A factor corresponding to the argument with the
    same name.}
  
}
\note{
  argument \code{minusMean} is now here as an "experimental"
  feature. The idea is that it could be used to detect
  non-stationarities of the reponses (in a repeated stimulation context)
  which would be correlated across different neurons. I'm not sure yet
  if this will be useful or not.
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com} }
\seealso{
  \code{\link{as.spikeTrain}},
  \code{\link{as.repeatedTrain}},
  \code{\link{data.frame}},
  \code{\link{factor}},
  \code{\link{rateEvolution}},
}
\examples{
## load Purkinje cell data recorded in cell-attached mode
data(sPK)
## coerce sPK to a spikeTrain object
sPK <- lapply(sPK, as.spikeTrain)
## get a rate evolution data frame
sPKreDF <- mkREdf(sPK)
## display result using coplot
coplot(rate ~ time | longitudinal,data=sPKreDF,panel=lines,show.given=FALSE)
\dontrun{
## make it prettier with with xyplot of package lattice
library(lattice)
xyplot(rate ~ time | longitudinal, data=sPKreDF,panel=panel.lines)
## if ggplot2 is installed, try it out
library(ggplot2)
qplot(time,rate,data=sPKreDF,geom="line",colour=longitudinal)
}

## load Purkinje cell data recorded with the NeuroNexus probes
data(mPK)
mPK <- lapply(mPK, as.repeatedTrain)
## get a rate evolution data frame
mPKreDF <- mkREdf(mPK)
## use coplot to display result
coplot(rate ~ time | longitudinal * across,data = mPKreDF,panel=lines)
\dontrun{
## make it prettier with with xyplot of package lattice
library(lattice)
xyplot(rate ~ time | across,data = mPKreDF,groups=longitudinal,panel=panel.lines)
xyplot(rate ~ time | across * longitudinal,data = mPKreDF, panel=panel.lines)
## if ggplot2 is installed, try it out
library(ggplot2)
qplot(time,rate,data=mPKreDF,geom="line",colour=longitudinal,facets=across ~ .)
}

## another example with the CAL1V data set
data(CAL1V)
CAL1V <- lapply(CAL1V,as.repeatedTrain)
## generate the data frame specifying the longitudinal argument
## to end up with a clearer display
CAL1VreDF <- mkREdf(CAL1V,longitudinal=paste(1:20))
coplot(rate ~ time | across * longitudinal,data=CAL1VreDF,panel=lines,show.given=FALSE)
\dontrun{
## if ggplot2 is installed, try it out
library(ggplot2)
qplot(time,rate,data=CAL1VreDF,geom="line",facets=longitudinal ~ across)
}

## another example with the CAL2C data set
data(CAL2C)
CAL2C <- lapply(CAL2C,as.repeatedTrain)
## generate the data frame specifying the longitudinal argument
## to end up with a clearer display
CAL2CreDF <- mkREdf(CAL2C,longitudinal=paste(1:20))
coplot(rate ~ time | across * longitudinal,data=CAL2CreDF,panel=lines,show.given=FALSE)
\dontrun{
## if ggplot2 is installed, try it out
library(ggplot2)
qplot(time,rate,data=CAL2CreDF,geom="line",facets=longitudinal ~ across)
}

}
\keyword{ts}
