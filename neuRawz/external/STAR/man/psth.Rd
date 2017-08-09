\name{psth}
\alias{psth}
\alias{plot.psth}
\title{Compute and Plot Peri-Stimulus Time Histogram}
\description{
  \code{psth} computes and \code{plot.psth} plots a peri-stimulus time histogram (called PST,
  post-stimulus time histogram by Gerstein and Kiang (1960)) from
  repeated presentations of a stimulation. Confidence bands can be
  obtained using the Poisson approximation.
}
\usage{
psth(repeatedTrain, breaks = 20, include.lowest = TRUE,
     right = TRUE, plot = TRUE, CI = 0.95, ...)
\method{plot}{psth}(x, stimTimeCourse = NULL, colStim = "grey80",
          colCI = NULL, xlab, ylab, main, xlim, ylim, lwd = 2,
          col = 1, ...)
}
\arguments{
  \item{repeatedTrain}{a \code{repeatedTrain} object or a list which can be
    coerced to such an object.}
  \item{x}{a \code{psth} object.}
  \item{stimTimeCourse}{\code{NULL} (default) or a two elements vector
    specifying the time boundaries (in s) of a stimulus presentation.}
  \item{colStim}{the background color used for the stimulus.}
  \item{breaks}{a numeric. A single number is interpreted has the number
    of bins; a vector of length 2 is interpreted as the bin width and
    the step to use (see details); otherwise interpreted as the
    position of the "breaks" between bins.}
  \item{include.lowest}{corresponding argument of \code{\link{hist}}.}
  \item{right}{corresponding argument of \code{\link{hist}}.}
  \item{plot}{corresponding argument of \code{\link{hist}}.}
  \item{CI}{The coverage probability of the confidence intervals.}
  \item{colCI}{if not \code{NULL} (default) a confidence band is
    plotted with the specified color; two dashed lines are plotted otherwise.}
  \item{xlim}{a numeric (default value supplied). See
    \code{\link{plot}}.}
  \item{ylim}{a numeric (default value supplied). See \code{\link{plot}}.}
  \item{xlab}{a character (default value supplied). See \code{\link{plot}}.}
  \item{ylab}{a character (default value supplied). See \code{\link{plot}}.}
  \item{main}{a character (default value supplied). See \code{\link{plot}}.}
  \item{lwd}{line width used to plot the estimated density. See \code{\link{plot}}.}
  \item{col}{color used to plot the estimated density. See \code{\link{plot}}.}
  \item{\dots}{see \code{\link{plot}}.}
}
\details{
  When confidence bands are requested they are obtained from the
  qunatiles of the \code{\link{Poisson}} distribution.

  When a 2 elements vector is used as \code{breaks} argument it is
  interpreted as specifying a bin width (first element if elements are
  unnamed, \code{"bw"} element otherwise) and a step (second element if
  elements are unnamed, \code{"step"} element otherwise). The idea is
  then to obtain a smoother looking PSTH by counting spikes within
  overlapping bins. That is if the center of the ith bin is xi the one
  of the (i+1)th bin will be xi + step.

}
\value{
  When \code{plot} is set to \code{FALSE} in \code{psth}, a list of
  class \code{psth} is returned and no plot
  is generated. This list has the following components:
  \item{freq }{a vector containing the instantaneous firing rate.}
  \item{ciUp }{a vector with the upper limit of the confidence band.}
  \item{ciLow }{a vector with the lower limit of the confidence band.}
  \item{breaks }{a numeric vector with the breaks in between which spikes
  were counted. Similar to the component of the same name returned by
  \code{\link{hist}}.}
  \item{mids }{a numeric vector with the mid points of
    \code{breaks}. Similar to the component of the same name returned
    by \code{\link{hist}}.}
  \item{counts }{a matrix with as many rows as components in
    \code{repeatedTrain} and as many columns as bins. Each element of
    the matrix contains the number of spikes falling in a given trial in
    a given bin.}
  \item{nbTrials }{the number of stimulations.}
  \item{call }{the matched call.}

  When \code{plot} is set to \code{TRUE} nothing is returned and a plot
  is generated as a side effect. Of course the same occurs upon calling
  \code{plot.psth} with a \code{psth} object argument.
  
}
\references{
  Gerstein, George L. and Kiang, Nelson Y.-S. (1960) An Approach to the
  Quantitative Analysis of Electrophysiological Data from Single
  Neurons. \emph{Biophysical Journal} \bold{1}:
  15--28. \url{http://www.pubmedcentral.nih.gov/articlerender.fcgi?tool=pubmed&pubmedid=13704760}

  Kalbfleisch, J. G. (1985) \emph{Probability and Statistical
    Inference. Volume 2: Statistical Inference}. Springer-Verlag.
}
\author{Christophe Pouzat \email{christophe.pouzat@gmail.com} }
\seealso{
  \code{\link{as.repeatedTrain}},
  \code{\link{is.repeatedTrain}},
  \code{\link{print.repeatedTrain}},
  \code{\link{summary.repeatedTrain}},
  \code{\link{raster}}
}
\examples{
## Load Vanillin responses data (first cockroach data set)
data(CAL1V)
## convert them into repeatedTrain objects
## The stimulus command is on between 4.49 s and 4.99s
CAL1V <- lapply(CAL1V,as.repeatedTrain)
## look at the individual raster plots
plot(CAL1V[["neuron 1"]],stimTimeCourse=c(4.49,4.99),main="N1")
## Create a simple black and white PSTH for neuron 1
psth(CAL1V[["neuron 1"]],stimTimeCourse=c(4.49,4.99),breaks=20)
## Rebuilt the same PSTH but with red confidence bands
psth(CAL1V[["neuron 1"]],stimTimeCourse=c(4.49,4.99),breaks=20,colCI=2)
## Make the PSTH smoother
psth(CAL1V[["neuron 1"]],stimTimeCourse=c(4.49,4.99),breaks=c(bw=0.5,step=0.05),colCI=2)
## Make a plot with PSTHs from 4 neurons superposed
## First get lists containing PSTHs from each neuron
psth1 <- psth(CAL1V[["neuron 1"]],breaks=c(bw=0.5,step=0.05),plot=FALSE)
psth2 <- psth(CAL1V[["neuron 2"]],breaks=c(bw=1,step=0.1),plot=FALSE)
psth3 <- psth(CAL1V[["neuron 3"]],breaks=c(bw=0.5,step=0.05),plot=FALSE)
psth4 <- psth(CAL1V[["neuron 4"]],breaks=c(bw=2,step=0.2),plot=FALSE)
## Get the maximal frequency to display
maxFreq <- max(max(psth1$ciUp),max(psth2$ciUp),max(psth3$ciUp),max(psth4$ciUp))
## Build plot
plot(c(0,10),c(0,75),type="n",
     xaxs="i",yaxs="i",xlab="Time (s)",
     ylab="Freq. (Hz)",
     main="PSTHs from 4 simultaneously recorded neurons",
     sub="20 stimulations with vanillin were used.")
## Add rectangle corresponding to stimulation command
rect(4.49,0,4.99,75,col="grey80",lty=0)
## Add the neurons PSTHs as confidence bands
polygon(c(psth1$mids,rev(psth1$mids)),c(psth1$ciLow,rev(psth1$ciUp)),col=1,border=NA)
polygon(c(psth2$mids,rev(psth2$mids)),c(psth2$ciLow,rev(psth2$ciUp)),col=2,border=NA)
polygon(c(psth3$mids,rev(psth3$mids)),c(psth3$ciLow,rev(psth3$ciUp)),col=3,border=NA)
polygon(c(psth4$mids,rev(psth4$mids)),c(psth4$ciLow,rev(psth4$ciUp)),col=4,border=NA)
legend(0.1,maxFreq,legend=paste("neuron",1:4),lty=1,col=1:4,bty="n")
}
\keyword{ts}
\keyword{survival}
