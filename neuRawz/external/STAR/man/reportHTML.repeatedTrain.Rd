\name{reportHTML.repeatedTrain}
\alias{reportHTML.repeatedTrain}
\title{Performs Basic Spike Train Analysis and Generates a Report in
  HTML Format from a repeatedTrain Object}
\description{
  Performs a "standard" analysis on a \code{repatedTrain} object, writes
  results to disk and generates a report in html format.
}
\usage{
\method{reportHTML}{repeatedTrain}(object, filename, extension = "html",
           directory = getwd(), Title, binSize = 0.025,
           method = c("gsspsth0", "gsspsth", "gampsth"),
           stimTimeCourse = NULL, colCI = 2,
           doTimeTransformation = TRUE, k = 100, bs = "tp",
           doGamCheck = FALSE, ...)
}
\arguments{
  \item{object}{a \code{repeatedTrain} object.}
  \item{filename}{a character string. The generic name of all the files
    (html, png as well as \code{R} data files which will be
    generated. See also \code{\link[R2HTML]{HTMLInitFile}}.}
  \item{extension}{see \code{\link[R2HTML]{HTMLInitFile}}.}
  \item{directory}{the full or relative path to the directory where the
    results are going to be stored. See also \code{\link[R2HTML]{HTMLInitFile}}.}
  \item{Title}{See \code{\link[R2HTML]{HTMLInitFile}}. If missing a default
    value baed on \code{filename} is provided.}
  \item{binSize}{See \code{\link{gsspsth}}, \code{\link{gsspsth0}},
    \code{\link{gampsth}}.}
  \item{method}{A character string, the name of the function used to
    generate the smooth psth, one of: \code{\link{gsspsth}}, \code{\link{gsspsth0}},
    \code{\link{gampsth}}.}
  \item{stimTimeCourse}{See \code{\link{plot.repeatedTrain}} and
    \code{\link{plot.gsspsth}}, \code{\link{plot.gsspsth0}}, \code{\link{plot.gampsth}}.}
  \item{colCI}{See \code{\link{plot.gsspsth}},
    \code{\link{plot.gsspsth0}}, \code{\link{plot.gampsth}}.}
  \item{doGamCheck}{Should function
    \code{\link[mgcv]{gam.check}} be used on the inhomogenous Poisson fit
    performed to obtain the smooth PSTH if \code{method} was set to \code{gampsth}?}
  \item{doTimeTransformation}{Should the estimated integrated intensity
    be used to perform a time transformation and generate Ogata's test
    plots?}
  \item{k, bs}{See \code{\link{gampsth}}.}
  \item{\dots}{Passed to \code{\link{gsspsth}}, \code{\link{gsspsth0}},
    \code{\link{gampsth}}.}
}
\details{
  A raster plot is added first to the report
  (\code{\link{plot.transformedTrain}}) with a smooth PSTH
  (\code{\link{gsspsth}}, \code{\link{gsspsth0}}, 
    \code{\link{gampsth}}.)
  superposed. The summary of the inhomogenous Poisson fit leading the
  smooth PSTH is added next together with a short summary describing
  how accurate the hypothesis of constant intensity/rate made during the
  pre-processing of the \code{repeatedTrain} was in view of the
  estimated rate. Check \code{\link{gsspsth}}, \code{\link{gsspsth0}},
    \code{\link{gampsth}} for details. A plot of the
  smooth PSTH with 95\% CI (approximate in the case of
  \code{\link{gampsth}}) is added.  If \code{doGamCheck} is
  set to \code{TRUE} and if \code{method} is set to \code{gampsth} a
  diagnostic plot for the  fitted inhomogenous
  Poisson model is added. If \code{doTimeTransformation} is set to
  \code{TRUE} the estimated integrated intensity is used to perform a
  time transformation and Ogata's test plots are generated.

  A \code{R} data file (\code{filename.rda}) is also generated with the
  following objects:
  \itemize{
    \item \code{PoissonF}: the \code{\link[gss]{gssanova}},
  \code{\link[gss]{gssanova0}} or \code{\link[mgcv]{gamObject}} object containing
    the result of the \code{\link[gss]{gssanova}},
  \code{\link[gss]{gssanova0}} or \code{\link[mgcv]{gam}} fit with the inhomogenous
    Poisson model.
    \item \code{Lambda}: the integrated intensity of \code{repeatedTrain}
    under the inhomogenous Poisson model hypothesis. If
    \code{doTimeTransformation} was set to \code{TRUE}.
    \item \code{fct}: the matched call.
  }
}
\value{
  Nothing is returned, an html file and figures in png format are
  written to disk together with the \code{R} variables generated during
  the analysis.
}
\author{Christophe Pouzat \email{christophe.pouzat@gmail.com} }
\seealso{
  \code{\link{as.repeatedTrain}},
  \code{\link{plot.repeatedTrain}},
  \code{\link{summary.repeatedTrain}},
  \code{\link{gsspsth}},
  \code{\link{gsspsth0}},
  \code{\link{gampsth}},
  \code{\link{transformedTrain}},
  \code{\link{plot.transformedTrain}},
  \code{\link{summary.transformedTrain}},
  \code{\link[gss]{gssanova}},
  \code{\link[gss]{gssanova0}},
  \code{\link[mgcv]{gam}},
  \code{\link[mgcv]{gam.check}},
  \code{\link{frt}}
}
\examples{
\dontrun{
## load e070528citronellal data set
data(e070528citronellal)
## make a standard analysis on the first neuron
reportHTML(e070528citronellal[["neuron 1"]],"e070528citronellalN1",stim=c(6.14,6.64))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{smooth}
\keyword{regression}
