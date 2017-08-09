\name{reportHTML.spikeTrain}
\alias{reportHTML.spikeTrain}
\title{Performs Basic Spike Train Analysis and Generates a Report in
  HTML Format from a spikeTrain Object}
\description{
  Performs a "standard" analysis on a \code{spikeTrain} object, computing some
  cross-correlation statistics if additional \code{spikeTrain} objects
  are provided, writes results to disk and generates a report in html format.
}
\usage{
\method{reportHTML}{spikeTrain}(object, filename, extension = "html",
           directory = getwd(), Title, forceTT = TRUE,
           digits = 3, timeUnit = "s", otherST,
           laglim = c(-0.1, 0.1),
           cch = c("both", "scch", "cch"),
           method = c("gsslockedTrain0", "gsslockedTrain", "gamlockedTrain"),
           doGamCheck = FALSE, k = 100, bs = "tp",
           nbEvtPerBin = 10, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a \code{spikeTrain} object.}
  \item{filename}{a character string. The generic name of all the files
    (html, png as well as \code{R} data files which will be
    generated. See also \code{\link[R2HTML]{HTMLInitFile}}.}
  \item{extension}{see \code{\link[R2HTML]{HTMLInitFile}}.}
  \item{directory}{the full or relative path to the directory where the
    results are going to be stored. See also \code{\link[R2HTML]{HTMLInitFile}}.}
  \item{Title}{See \code{\link[R2HTML]{HTMLInitFile}}. If missing a default
    value baed on \code{filename} is provided.}
  \item{forceTT}{Should a time transformation be performed and the
    \code{\link{compModels}} plots be generated even if none of the six
    renewal models fits the data?}
  \item{timeUnit, digits}{see \code{\link{summary.spikeTrain}}.}
  \item{otherST}{a named list of \code{spikeTrain} objects from
    simultaneously recorded neurons or nothing.}
  \item{laglim}{see \code{\link{lockedTrain}}.}
  \item{cch}{if \code{otherST} is given (ie, not missing)
    cross-intensity plots will be made using the neuron of
    \code{spikeTrain} as a reference. Should smooth version of the
    cross-intensity be computed (\code{"scch"}), a "classical" one
    (\code{"cch"}) or both (\code{"both"}). Only the first element of
    \code{cch} is used.}
  \item{method}{A character string, the name of the function used to
    generate the smooth cross-correlation histograms, one of:
    \code{\link{gsslockedTrain}}, 
    \code{\link{gsslockedTrain0}}, 
    \code{\link{gamlockedTrain}}.}
  \item{doGamCheck}{if smooth estimates are requested and \code{method}
    is set to \code{\link{gamlockedTrain}}, should function
    \code{\link[mgcv]{gam.check}} be used on them?}
  \item{k}{see \code{\link{gamlockedTrain}}.}
  \item{bs}{see \code{\link{gamlockedTrain}}.}
  \item{nbEvtPerBin}{a number of event per bin used in a way similar to
    the argument with the same name in \code{\link{jpsth}} when a bining
  is used for pre-processing.}
  \item{\dots}{Passed to \code{\link{gsslockedTrain}},
    \code{\link{gsslockedTrain0}}, 
    \code{\link{gamlockedTrain}}.}
}
\details{
  A spike train plot (\code{\link{plot.spikeTrain}}) is performed
  first. The summary (\code{\link{summary.spikeTrain}}) is computed next
  and part of its output is written to the html file. The renewal tests
  are then carried out and their results added
  (\code{\link{renewalTestPlot}}). The six duration distributions are
  fitted (\code{\link{compModels}} with argument \code{plot} set to
  \code{FALSE}) and the best one is used to apply a time transformation
  to \code{spikeTrain}. The Ogata's tests are applied
  (\code{\link{summary.transformedTrain}}) and if they are all within
  the 99\% confidence interval, the result of the transformation is
  plotted (\code{\link{plot.transformedTrain}}) as well as all the Q-Q
  plots of \code{\link{compModels}}. If \code{forceTT} is set to
  \code{TRUE} (default), then these last two plots are added even if the
  best model does not pass the tests.

  If other \code{spikeTrain} objects are provided as a named list via
  argument \code{otherST}, then cross-correlation/cross-intensity
  functions are estimated; Two estimations methods are available, the
  classical histogram and a smooth version of it. Argument \code{cch}
  controls if a single estimation is performed or if both are
  performed. If the smooth version is requested a summary of the
  \code{\link[gss]{gssanova}}, \code{\link[gss]{gssanova0}} or
  \code{\link[mgcv]{gam}} fit is printed (depending on the chosen value
  for argument \code{method}). Moreover if argument
  \code{doGamCheck} is set to \code{TRUE} (and if \code{method} is set
  to \code{\link{gamlockedTrain}}) then check plots
  (\code{\link[mgcv]{gam.check}}) are added to the report.

  A \code{R} data file (\code{filename.rda}) is also generated with the
  following objects:
  \itemize{
    \item \code{cm}: the result of \code{\link{compModels}}.
    \item \code{bestFit}: the \code{durationDistribution} object
    returned obtained by fitting the best model among the 6.
    \item \code{Lambda}: the integrated intensity of \code{spikeTrain}
    with the best model.
    \item \code{fct}: the matched call.
    \item \code{cchL}: if other trains were provided and if argument
    \code{cch} was set to \code{"both"} or to \code{"cch"}. A list with
    as many components as the \code{otherST} argument. Each component is
    the a \code{hist.lockedTrain} object.
    \item \code{scchL}: if other trains were provided and if argument
    \code{cch} was set to \code{"both"} or to \code{"scch"}. A list with
    as many components as the \code{otherST} argument. Each component is
    the a \code{gsslockedTrain}, \code{gsslockedTrain0} or
    \code{gamlockedTrain} object.
  }
  
}
\value{
  Nothing is returned, an html file and figures in png format are
  written to disk together with the \code{R} variables generated during
  the analysis.
}
\author{Christophe Pouzat \email{christophe.pouzat@gmail.com} }
\seealso{
  \code{\link{as.spikeTrain}},
  \code{\link{plot.spikeTrain}},
  \code{\link{summary.spikeTrain}},
  \code{\link{renewalTestPlot}},
  \code{\link{plot.spikeTrain}},
  \code{\link{compModels}},
  \code{\link{transformedTrain}},
  \code{\link{plot.transformedTrain}},
  \code{\link{summary.transformedTrain}},
  \code{\link[gss]{gssanova}},
  \code{\link[gss]{gssanova0}},
  \code{\link[mgcv]{gam}},
  \code{\link[mgcv]{gam.check}},
  \code{\link{lockedTrain}},
  \code{\link{gsslockedTrain}},
  \code{\link{gsslockedTrain0}},
  \code{\link{gamlockedTrain}}
}
\examples{
\dontrun{
## load e070528spont data set
data(e070528spont)
## perform a standard analysis on neuron 1, looking for cross-correlations
## with the 3 other neurons up to lag +/- 250 ms.
## Store the results under the generic name: e070528spontN1
reportHTML(e070528spont[["neuron 1"]],"e070528spontN1",otherST=e070528spont[-1],laglim=c(-1,1)*0.25,forceTT=FALSE)
## Neuron 1 of e070528spont is exceptional in that it can be well
## described by a renewal process...}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{smooth}
\keyword{regression}
