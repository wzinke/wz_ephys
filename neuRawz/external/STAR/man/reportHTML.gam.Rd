\name{reportHTML.gam}
\alias{reportHTML.gam}
\title{Generates a Report in
  HTML Format from a STAR gam Object }
\description{
  Writes the result of a \code{gam} fit in an html file.
}
\usage{
\method{reportHTML}{gam}(object, filename, extension = "html",
           directory = getwd(), Title,
           neuron, neuronEvts, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object returned by \code{\link[mgcv]{gam}}.}
  \item{filename}{a character string. The generic name of all the files
    (html, png as well as \code{R} data files which will be
    generated. See also \code{\link[R2HTML]{HTMLInitFile}}.}
  \item{extension}{see \code{\link[R2HTML]{HTMLInitFile}}.}
  \item{directory}{the full or relative path to the directory where the
    results are going to be stored. See also \code{\link[R2HTML]{HTMLInitFile}}.}
  \item{Title}{See \code{\link[R2HTML]{HTMLInitFile}}. If missing a default
    value baed on \code{filename} is provided.}
  \item{neuron}{a character string describing to which the analysis
    refers and used for the titles of the interaction plots (see \code{\link{plot.frt}}).}
  \item{neuronEvts}{a named list with the \code{event} variable from the
  data frame returned by \code{\link{mkGLMdf}} and corresponding to the
  other neurons recorded simultaneously. One list element per neuron.}
  \item{\dots}{Not used, only there for compatibilty with the generic
    method definition.}
}
\details{
  A summary (\code{\link[mgcv]{summary.gam}}) of \code{object} is added
  to the report. A plot of the spike train after time transformation
  \code{\link{transformedTrain}} comes next followed by a renewal test
  plot (\code{\link{renewalTestPlot}}) of the spike train on the time
  transformed scale. The "usual" Ogata's tests plots
  (\code{\link{plot.transformedTrain}}) are added. Then if other trains
  are provided as a named list via argument \code{neuronEvts},
  interactions plots (\code{\link{plot.frt}}) are built showing both the
  survivor function and the Berman's test. The report ends with the
  \code{call} which generated \code{object}. 
}
\value{
  Nothing is returned, an html file and figures in png format are
  written to disk.
}
\author{Christophe Pouzat \email{christophe.pouzat@gmail.com} }
\seealso{
  \code{\link{mkGLMdf}},
  \code{\link[mgcv]{gam}},
  \code{\link[mgcv]{gam.check}},
  \code{\link{frt}},
  \code{\link{transformedTrain}},
  \code{\link{plot.transformedTrain}},
  \code{\link{summary.transformedTrain}}
}
\examples{
\dontrun{
## load e070528spont data set
data(e070528spont)
## make a data frame for gam using a 2 ms bin width 
spontDF <- mkGLMdf(e070528spont,0.002,0,60)
## make data frames specific of each neuron
n1.spontDF <- spontDF[spontDF$neuron=="1",]
n2.spontDF <- spontDF[spontDF$neuron=="2",]
n3.spontDF <- spontDF[spontDF$neuron=="3",]
n4.spontDF <- spontDF[spontDF$neuron=="4",]
## save space by removing the now redundant spontDF
rm(spontDF)
## fit neuron 1 using the gam representation of a
## renewal process and a binomial model
n1.spontFit1 <- gam(event ~ s(lN.1,k=25,bs="cr"),data=n1.spontDF,family=binomial())
## create a list with the discretized spike times of the 3 other neurons
preN1 <- list(n2=with(n2.spontDF,event),n3=with(n3.spontDF,event),n4=with(n4.spontDF,event))
## generate the report
reportHTML(n1.spontFit1,"e070528spontN1gFit",neuron="1",neuronEvts=preN1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{smooth}
\keyword{regression}
