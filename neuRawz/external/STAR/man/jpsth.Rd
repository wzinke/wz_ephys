\name{jpsth}
\alias{jpsth}
\alias{jsd}
\alias{contour.jpsth}
\alias{image.jpsth}
\alias{persp.jpsth}
\alias{jpsth2df}

\title{Related Functions and Methods for Joint-PSTHs and Joint Scatter Diagrams}
\description{
  Some mainly graphical tools to probe interactions between 2 neurons
  recorded in the presence of a repeated stimulation.
}
\usage{
jsd(xRT, yRT, acquisitionWindow, xlab, ylab,
    main, pch = ".", ...)
jpsth(xRT, yRT, xBreaks, yBreaks,
      acquisitionWindow, nbEvtPerBin = 50)
\method{contour}{jpsth}(x, xlab, ylab, main, ...)
\method{image}{jpsth}(x, xlab, ylab, main, ...)
\method{persp}{jpsth}(x, xlab, ylab, main, ...)
jpsth2df(object)
}

\arguments{
  \item{xRT}{a \code{repeatedTrain} object whose spike times will
  appear on the abscissa of the plots.}
  \item{yRT}{a \code{repeatedTrain} object whose spike times will
  appear on the ordinate of the plots. It must have the same length as \code{xRT}.}
  \item{x, object}{\code{jpsth} objects.}
  \item{xBreaks, yBreaks}{A single number (the bin width) or a vector defining
    bins boundaries on the X and Y axis. If missing a default is provided.}
  \item{acquisitionWindow}{a 2 elements vector specifying the begining
    and the end of the acquisition. If missing values are obtained using
    the  \code{\link{floor}} of the smallest spike time and the
    \code{\link{ceiling}} of the largest one.}
  \item{nbEvtPerBin}{If both \code{xBreaks} and \code{xBreaks} are
    missing a bin width, \code{bw}, is computed such that the expected
    value of the count per cell (2 dimensional bin) would be
    \code{nbEvtPerBin} assuming a stationary Poisson discharge for both
    neurons.}
  \item{xlab}{a character (default value supplied). See \code{\link{plot}}.}
  \item{ylab}{a character (default value supplied). See \code{\link{plot}}.}
  \item{main}{a character (default value supplied). See \code{\link{plot}}.}
  \item{pch}{the type of "points" displayed by \code{jsd}. See \code{\link{plot}}.}
  \item{\dots}{additional arguments passed to \code{\link{plot}} by
    \code{jsd} and to respective generic methods by
    \code{contour.jpsth}, \code{image.jpsth} and \code{persp.jpsth}.}
}
\details{
  The joint scatter diagram was introduced by Gerstein and Perkel
  (1972). The joint peristimulus time histogram is a binned version of
  it (Aertsen et al, 1989). \code{jpsth2df} allows the reformating of a \code{jpsth} object in
  order to compute a smooth version of it with
  \code{\link[gss]{gssanova}}, \code{\link[gss]{gssanova0}} or \code{\link[mgcv]{gam}}.
}
\value{
  \code{jsd} is used for its side effect, a plot is generated and
  nothing is returned.

  \code{jpsth2df} returns a \code{\link{data.frame}} with the following
  variables: \code{Count}, the counts per cell; \code{X}, the position
  of the cell on the X axis; \code{Y}, the position of the cell on the Y
  axis; 
  and attributes: \code{xBreaks}, \code{yBreaks}, \code{xTotal},
  \code{yTotal}, \code{nbTrials}, \code{acquisitionWindow} corresponding
  to the components of its argument with the same name and
  \code{originalCall} corresponding to component \code{call}.
  
  \code{jpsth} returns a list of class \code{jpsth} with the following
  components: 
  \item{counts }{a matrix storing the counts per cell.}
  \item{density }{a matrix storing the density in each cell.}
  \item{xMids }{a vector containing the X positions of the cells.}
  \item{yMids }{a vector containing the Y positions of the cells.}
  \item{xBreaks }{a vector containing the bin boundaries of the cells
    along the X axis.}
  \item{yBreaks }{a vector containing the bin boundaries of the cells
    along the X axis.}
  \item{xTotal }{the total number of spikes of the "X" neuron.}
  \item{yTotal }{the total number of spikes of the "Y" neuron.}
  \item{xFreq }{the mean freqency of the "X" neuron.}
  \item{yFreq }{the mean freqency of the "Y" neuron.}
  \item{nbTrials }{the number of trials of \code{xRT} (and \code{yRT}).}
  \item{acquisitionWindow }{the boundaries of the acquisition window.}
  \item{call }{the matched call.}

  
}
\references{
  Gerstein, G. L. and Perkel, D. H. (1972) Mutual temporal relationships
  among neuronal spike trains. Statistical techniques for display and
  analysis. \emph{Biophys J} \bold{12}:
  453--473. \url{http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=1484144} 

  Aertsen, A. M., Gerstein, G. L., Habib, M. K., Palm, G. (1989)
  Dynamics of neuronal firing correlation: modulation of "effective
  connectivity". \emph{J Neurophysiol} \bold{61}:
  900--917. \url{http://jn.physiology.org/cgi/content/abstract/61/5/900} 
}
\author{Christophe Pouzat \email{christophe.pouzat@gmail.com} }
\note{
  I use "joint scatter diagram" for what Gerstein and Perkel (1972) more
  properly call a "joint peristimulus time scatter diagram".
}
\seealso{
  \code{\link{lockedTrain}},
  \code{\link{plot.lockedTrain}},
  \code{\link{hist.lockedTrain}},
  \code{\link{gsslockedTrain}},
  \code{\link{plot.gsslockedTrain}},
  \code{\link{gsslockedTrain0}},
  \code{\link{plot.gsslockedTrain0}},
  \code{\link{gamlockedTrain}},
  \code{\link{plot.gamlockedTrain}},
  \code{\link{contour}},
  \code{\link{image}},
  \code{\link{persp}},
  \code{\link{attr}},
  \code{\link{attributes}}
}
\examples{
## load e070528citronellal data
data(e070528citronellal)
## plot a jsd with neuron 1 on X and neuron 2 on Y
jsd(e070528citronellal[[1]],e070528citronellal[[2]])
## now make the jpsth
j1.2 <- jpsth(e070528citronellal[[1]],e070528citronellal[[2]])
## make a contour plot
contour(j1.2)
## make an image plot
image(j1.2)
## make a persp plot
persp(j1.2)
\dontrun{
## fit a gss model with interactions
## use a larger bin width for the jpsth
j1.2 <- jpsth(e070528citronellal[[1]],e070528citronellal[[2]],0.2,0.2)
## get a data frame
j1.2DF <- jpsth2df(j1.2)
## To save computation time start analyzing
## just before the stimulation time
j1.2DF <- j1.2DF[j1.2DF$X > 6 & j1.2DF$Y>6,]
gf <- gssanova(Count ~ X*Y, family="poisson", data=j1.2DF,seed=20061001)
## Use the project function of gss to check the significance
## of the interaction term
project(gf2,inc=c("X","Y"))
}
\dontrun{
## fit a gam model assuming no interaction
## get a data frame
j1.2DF <- jpsth2df(j1.2)
fitNoI <- gam(Count ~ s(X,k=100,bs="cr") + s(Y,k=100,bs="cr"),data=j1.2DF,family=poisson())
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
