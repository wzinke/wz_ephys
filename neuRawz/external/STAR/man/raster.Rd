\name{raster}
\alias{raster}
\alias{plot.repeatedTrain}
\title{ Generate a Raster Plot }
\description{
  Given a list of spike trains (or a \code{repeatedTrain} object) where
  each train was acquired during,
  say, one presentation of a given stimulus, a raster plot is
  generated. If stimulus time properties are specified, the stimulus
  application time also appears on the plot. 
}
\usage{
\method{plot}{repeatedTrain}(x, stimTimeCourse = NULL,
       colStim = "grey80", xlim, pch, xlab, ylab, main, ...)
raster(x, stimTimeCourse = NULL, colStim = "grey80",
       xlim, pch, xlab, ylab, main, ...)
}
\arguments{
  \item{x}{a \code{repeatedTrain} object or a list which can be
    coerced to such an object.}
  \item{stimTimeCourse}{\code{NULL} (default) or a two elements vector
    specifying the time boundaries (in s) of a stimulus presentation.}
  \item{colStim}{the background color used for the stimulus.}
  \item{xlim}{a numeric (default value supplied). See \code{\link{plot}}.}
  \item{pch}{data symbol used for the spikes. See \code{\link{plot}}.}
  \item{xlab}{a character (default value supplied). See \code{\link{plot}}.}
  \item{ylab}{a character (default value supplied). See \code{\link{plot}}.}
  \item{main}{a character (default value supplied). See \code{\link{plot}}.}
  \item{\dots}{see \code{\link{plot}}.}
}
\details{
  Basic raster plot stuff.
}
\value{
  Nothing is returned \code{raster} is used for its side effect, a
  plot is generated on the current graphical device.
}
\references{
  Brillinger, David R. (1992) Nerve Cell Spike Train Data
  Analysis: A Progression of Technique. \emph{JASA} \bold{87}: 260--271. 
}
\author{Christophe Pouzat \email{christophe.pouzat@gmail.com} }
\note{
  Brillinger (1992) calls these plots "rastor" instead of raster...
}
\seealso{
  \code{\link{as.repeatedTrain}},
  \code{\link{is.repeatedTrain}},
  \code{\link{print.repeatedTrain}},
  \code{\link{summary.repeatedTrain}},
  \code{\link{psth}}
}
\examples{
## Load Vanillin responses data (first cockroach data set)
data(CAL1V)
## convert them into repeatedTrain objects
## The stimulus command is on between 4.49 s and 4.99s
CAL1V <- lapply(CAL1V,as.repeatedTrain)
## look at the individual raster plots
raster(CAL1V[["neuron 1"]],stimTimeCourse=c(4.49,4.99),main="N1")
plot(CAL1V[["neuron 2"]],stimTimeCourse=c(4.49,4.99),main="N2")
plot(CAL1V[["neuron 3"]],stimTimeCourse=c(4.49,4.99),main="N3")
plot(CAL1V[["neuron 4"]],stimTimeCourse=c(4.49,4.99),main="N4")
}
\keyword{ts}
\keyword{survival}
