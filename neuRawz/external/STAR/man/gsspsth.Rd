\name{gsspsth}
\alias{gsspsth}
\alias{gsspsth0}
\alias{print.gsspsth}
\alias{print.gsspsth0}
\alias{summary.gsspsth}
\alias{summary.gsspsth0}
\alias{plot.gsspsth}
\alias{plot.gsspsth0}
\alias{simulate.gsspsth0}
\alias{simulate.gsspsth}
\title{Smooth Peri Stimulus Time Histogram Related Functions and
  Methods: The Smoothing Spline Approach}
\description{
  Function \code{gsspsth} and \code{gsspsth0} compute a smooth psth, while method
  \code{print.gsspsth} and \code{print.gsspsth0} print and
  \code{summary.gsspsth} or \code{summary.gsspsth0} summarize the
  \code{gssanova} /  \code{gssanova0} objects contained in the returned \code{gsspsth} or
  \code{gsspsth0} objects,
  \code{plot.gsspsth} or \code{plot.gsspsth0} plot them and
  \code{simulate.gsspsth} or \code{simulate.gsspsth0} simulate data from
  fitted objects. 
}
\usage{
gsspsth(repeatedTrain, binSize = 0.025, plot = FALSE, ...)
gsspsth0(repeatedTrain, binSize = 0.025, plot = FALSE, ...)
\method{print}{gsspsth}(x, ...)
\method{print}{gsspsth0}(x, ...)
\method{summary}{gsspsth}(object, ...)
\method{summary}{gsspsth0}(object, ...)
\method{plot}{gsspsth}(x, stimTimeCourse = NULL, colStim = "grey80",
           colCI = NULL, xlab, ylab, main, xlim, ylim,
           lwd = 2, col = 1, ...)
\method{plot}{gsspsth0}(x, stimTimeCourse = NULL, colStim = "grey80",
           colCI = NULL, xlab, ylab, main, xlim, ylim,
           lwd = 2, col = 1, ...)
\method{simulate}{gsspsth}(object, nsim = 1, seed = NULL, ...)
\method{simulate}{gsspsth0}(object, nsim = 1, seed = NULL, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{repeatedTrain}{a \code{repeatedTrain} object or a list which can be
    coerced to such an object.}
  \item{binSize}{the bin size (in s) used to generate the observations on which
  the gss fit will be performed. See details below.}
  \item{plot}{corresponding argument of \code{\link{hist}}. Should a
    plot be generated or not?}
  \item{object}{a \code{gsspsth} or a \code{gsspsth0} object.}
  \item{x}{a \code{gsspsth} or a \code{gsspsth0} object.}
  \item{stimTimeCourse}{\code{NULL} (default) or a two elements vector
    specifying the time boundaries (in s) of a stimulus presentation.}
  \item{colStim}{the background color used for the stimulus.}
  \item{colCI}{if not \code{NULL} (default) a confidence band is
    plotted with the specified color; two dashed lines are plotted otherwise.}
  \item{xlim}{a numeric (default value supplied). See
    \code{\link{plot}}.}
  \item{ylim}{a numeric (default value supplied). See \code{\link{plot}}.}
  \item{xlab}{a character (default value supplied). See \code{\link{plot}}.}
  \item{ylab}{a character (default value supplied). See \code{\link{plot}}.}
  \item{main}{a character (default value supplied). See \code{\link{plot}}.}
  \item{lwd}{line width used to plot the estimated density. See \code{\link{plot}}.}
  \item{col}{color used to plot the estimated density. See
    \code{\link{plot}}.}
  \item{nsim}{number of \code{repeatedTrain} objects to simulate. Defaults to 1.}
  \item{seed}{see \code{\link{simulate}}.}
  \item{\dots}{in \code{gsspsth}, respectively \code{gsspsth0}, the
    \dots are passed to the internally called \code{\link[gss]{gssanova}}, repectively
    \code{\link[gss]{gssanova0}}. In
    \code{plot.gsspsth} and \code{plot.gsspsth0} they are passed to
    \code{\link{plot}} which is 
    called internally. They are not used otherwise.}
}
\details{
  \code{gsspsth} calls internally \code{\link[gss]{gssanova}} while
  \code{gsspsth0} calls \code{\link[gss]{gssanova0}}. See the respective
  documentations and references therein for an explanation of the differences.
  For both \code{gsspsth} and \code{gsspsth0}, the raw data contained in
  \code{repeatedTrain} are 
  pre-processed with \code{\link{hist}} using a bin size given by
  argument \code{binSize}. This \code{binSize} should be small "enough". That is, the
  rate of the aggregated train created by collapsing the spike times of
  the different trials onto a single "pseudo" spike train, should not
  change too much on the scale of \code{binSize} (see Ventura et al
  (2002) Sec. 4.2 p8 for more details). Argument \code{nbasis} of
  \code{\link[gss]{gssanova}} called internally by \code{gsspsth} is set
  to the number of bins of the histogram resulting from the
  preprocessing stage.

  \code{simulate.gsspsth} and \code{simulate.gsspsth0} perform exact
  simuations of inhomogenous Poisson processes whose (time dependent)
  rate/intensity function is accessible through the componenent
  \code{lambdaFct} of objects of class \code{gsspsth} and
  \code{gsspsth0}. The inhomogenous Poisson processes are simulated with
  the thinning method (Devroye, 1986, pp 253-256).
}
\value{
  When \code{plot} is set to \code{FALSE} in \code{gsspsth}, repectively
  \code{gsspsth0}, a list of
  class \code{gsspsth}, respectively \code{gsspsth0}, is returned and no plot
  is generated. These list have the following components:
  \item{freq }{a vector containing the instantaneous firing rate in the
    middle of the "thin" bins used for preprocessing.}
  \item{ciUp }{a vector with the upper limit of a pointwise 95\% confidence interval.   Check \code{\link[gss]{predict.ssanova}} for details.}  
  \item{ciLow }{a vector with the lower limit of a pointwise 95\% confidence interval.} 
  \item{breaks }{a vector with 2 elements the ealiest and the latest spike in \code{repeatedTrain}.}
  \item{mids }{a numeric vector with the mid points of
    the bins.}
  \item{counts }{a vector with the actual number of spikes in each bin.}
  \item{nbTrials }{the number of trials in \code{repeatedTrain}.}
  \item{lambdaFct }{a function of a single time argument returning the
    estimated intensity (or instantaneous rate) at its argument.}
  \item{LambdaFct }{a function of a single time argument returning the
    integrale of estimated intensity (or instantaneous rate) at its
    argument. That is, the integrated intensity. \code{\link{integrate}}
  is used by this function.}
  \item{call }{the matched call.}

  When \code{plot} is set to \code{TRUE} nothing is returned and a plot
  is generated as a side effect. Of course the same occurs upon calling
  \code{plot.gsspsth} with a \code{gsspsth} object argument or
  \code{plot.gsspsth0} with a \code{gsspsth0}. 

  \code{print.gsspsth} returns the result of \code{\link[gss]{print}}
  applied to the \code{gssanova} object generated by \code{gsspsth}
  and stored in the \code{\link{environment}} of both \code{lambdaFct}
  and \code{LambdaFct}. The same goes for \code{print.gsspsth0}.

  \code{summary.gsspsth} returns the result of \code{\link[gss]{summary.gssanova}}
  applied to the \code{gssanova} object generated by \code{gsspsth}
  and stored in the \code{\link{environment}} of both \code{lambdaFct}
  and \code{LambdaFct}. The same goes for \code{summary.gsspsth0}.
  
  \code{simulate.gsspsth} and \code{simulate.gsspsth0} return a
  \code{repeatedTrain} object if argument \code{nsim} is set to one and
  a list of such objects if it is greater than one.
  
}
\references{
  Gu C. (2002) \emph{Smoothing Spline ANOVA Models}. Springer.
  
  Ventura, V., Carta, R., Kass, R. E., Gettner, S. N. and Olson,
  C. R. (2002) Statistical analysis of temporal evolution in
  single-neuron firing rates. \emph{Biostatistics} \bold{3}: 1--20.

  Kass, R. E., Ventura, V. and Cai, C. (2003) Statistical smoothing of
  neuronal data. \emph{Network: Computation in Neural Systems}
  \bold{14}: 5--15.

  Devroye Luc (1986) \emph{Non-Uniform Random Variate
    Generation}. Springer. Book available in pdf format at:
    \url{http://cg.scs.carleton.ca/~luc/rnbookindex.html}. 
}
\author{Christophe Pouzat \email{christophe.pouzat@gmail.com} }
\note{
  Most of the components of the list returned by \code{gsspsth} and
    \code{gsspsth0} are not of
  direct interest for the user but they are used by, for instance,
  \code{\link{reportHTML.repeatedTrain}}. 
}
\seealso{
  \code{\link{psth}},
  \code{\link{plot.psth}},
  \code{\link[gss]{gssanova}},
  \code{\link[gss]{gssanova0}},
  \code{\link[gss]{summary.gssanova}},
  \code{\link[gss]{summary.gssanova0}},
  \code{\link{reportHTML.repeatedTrain}},
  \code{\link{simulate}}
}
\examples{
\dontrun{
## Get the e070528citronellal data set into workspace
data(e070528citronellal)
## Compute gsspsth without a plot for neuron 1
## using a smmothing spline with gssanova0, and default bin size of 25 ms.
n1CitrGSSPSTH0 <- gsspsth0(e070528citronellal[[1]])
## plot the result
plot(n1CitrGSSPSTH0,stim=c(6.14,6.64),colCI=2)
## get a summary of the gss fit
summary(n1CitrGSSPSTH0)
## Now take a look at the observation on which the gss
## was actually performed
plot(n1CitrGSSPSTH0$mids,n1CitrGSSPSTH0$counts,type="l")
## Add the estimated smooth psth after proper scaling
theBS <- diff(n1CitrGSSPSTH0[["mids"]])[1]
Y <- n1CitrGSSPSTH0$lambdaFct(n1CitrGSSPSTH0$mids)*theBS*n1CitrGSSPSTH0$nbTrials
lines(n1CitrGSSPSTH0$mids,Y,col=4,lwd=2)

## check the (absence of) effect of the pre-binning by using a smaller
## and larger one, say 10 and 75 ms
n1CitrGSSPSTH0_10 <- gsspsth0(e070528citronellal[[1]],binSize=0.01)
n1CitrGSSPSTH0_75 <- gsspsth0(e070528citronellal[[1]],binSize=0.075)
## plot the "high resolution" smoothed-psth
plot(n1CitrGSSPSTH0_10,colCI="grey50")
## add to it the estimate obtained with the "low resolution" one
Y_75 <- n1CitrGSSPSTH0_75$lambdaFct(n1CitrGSSPSTH0_10$mids)
lines(n1CitrGSSPSTH0_10$mids,Y_75,col=2,lwd=2)

## simulate data from the first fitted model
s1 <- simulate(n1CitrGSSPSTH0)
## look at the result
s1

## Do the same thing with gsspsth
n1CitrGSSPSTH <- gsspsth(e070528citronellal[[1]])
plot(n1CitrGSSPSTH,stim=c(6.14,6.64),colCI=2)
summary(n1CitrGSSPSTH)
plot(n1CitrGSSPSTH$mids,n1CitrGSSPSTH$counts,type="l")
theBS <- diff(n1CitrGSSPSTH[["mids"]])[1]
Y <- n1CitrGSSPSTH$lambdaFct(n1CitrGSSPSTH$mids)*theBS*n1CitrGSSPSTH$nbTrials
lines(n1CitrGSSPSTH$mids,Y,col=4,lwd=2)
## check the (absence of) effect of the pre-binning by using a smaller
## and larger one, say 10 and 75 ms
n1CitrGSSPSTH_10 <- gsspsth(e070528citronellal[[1]],binSize=0.01)
n1CitrGSSPSTH_75 <- gsspsth(e070528citronellal[[1]],binSize=0.075)
## plot the "high resolution" smoothed-psth
plot(n1CitrGSSPSTH_10,colCI="grey50")
## add to it the estimate obtained with the "low resolution" one
Y_75 <- n1CitrGSSPSTH_75$lambdaFct(n1CitrGSSPSTH_10$mids)
lines(n1CitrGSSPSTH_10$mids,Y_75,col=2,lwd=2)
## simulate data from the first fitted model
s1 <- simulate(n1CitrGSSPSTH)
## look at the result
s1
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{smooth}
\keyword{regression}
