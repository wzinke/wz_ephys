\name{plot.quickPredict}
\Rdversion{1.1}
\alias{plot.quickPredict}
\alias{lines.quickPredict}
\alias{image.quickPredict}
\alias{contour.quickPredict}
\alias{persp.quickPredict}

\title{
Graphical Methods for quickPredict Objects
}
\description{
\code{plot}, \code{lines}, \code{image}, \code{contour} and \code{persp} methods for
\code{\link{quickPredict}} objects.
}
\usage{
\method{plot}{quickPredict}(x, y, xFct = function(x) x, style = c("band",
     "simple"), ylim, meanCol = 2, meanWD = 2, bandCol = "grey50", xlab,
     ylab, ...)
\method{lines}{quickPredict}(x, what = c("mean", "upr", "lwr"),
      xFct = function(x) x, ...)
\method{image}{quickPredict}(x, main, xlab, ylab, ...)
\method{contour}{quickPredict}(x, what = c("mean", "sd"), main, xlab, ylab, add,
        ...)
\method{persp}{quickPredict}(x, what = c("mean", "sd"), main, xlab, ylab, zlab, ...)
}
\arguments{
  \item{x}{ a \code{\link{quickPredict}} object.}
  \item{y}{ Not used but required by the plot method definition.}
  \item{xFct}{ a function applied to the \code{xx} component of argument
    \code{x} which is itself a \code{\link{quickPredict}} object. Useful
    to get a plot with the "native" scale of the variable.}
  \item{style}{ a character string, either \code{"band"} or \code{"simple"}. Controls the way
    confidence intervals are displayed. If \code{"band"} is selected
    they are shown as a ribbon. Their boundaries appear as
    dashed lines otherwise.}
  \item{ylim}{ see \code{\link{plot}} (default provided if missing).}
  \item{meanCol}{ the color used to display the estimated mean of the term
    (see argument \code{col} in \code{\link{plot}}).}
  \item{meanWD}{ the width used to display the estimated mean of the term
    (see argument \code{lwd} in \code{\link{plot}}).}
  \item{bandCol}{ the color of the confidence interval ribbon when one is
    is drawn.}
  \item{xlab, ylab, zlab}{ see \code{\link{plot}} and
    \code{\link{persp}} (default provided if any is missing).}
 \item{main}{ see \code{\link{image}}, \code{\link{contour}} and
   \code{\link{persp}} (default provided if missing).}
 \item{what}{ for \code{lines.quickPredict}, one of the following character strings: \code{"mean"}, \code{"upr"}, \code{"lwr"}.
   Controls the line drawn: the estimated mean, upper bound of
   the 95\% CI or lower bound. For \code{contour.quickPredict} and
   \code{persp.quickPredict}, a character string specifying if the mean
   or the sd contours or surface
   should be drawn. }
 \item{add}{ see \code{\link{contour}}.}
 \item{\dots}{ additional arguments to \code{\link{plot}}, \code{\link{lines}},
   \code{\link{image}}, \code{\link{contour}} and \code{\link{persp}}.}
}
\value{
Nothing returned, these functions are used for their side effects:
plots are generated or modified.
}
\author{
Christophe Pouzat \email{christophe.pouzat@gmail.com}
}

\seealso{
  \code{\link{quickPredict}},
  \code{\link{plot}},
  \code{\link{lines}}
  \code{\link{image}},
  \code{\link{contour}},
  \code{\link{persp}},
  \code{\link{plot.ssanova}}
}
\examples{
\dontrun{
## Follow up of ssanova example of gss
data(nox)
nox.fit <- ssanova(log10(nox)~comp*equi,data=nox)
## get prediction for the first term, comp
comp.pred <- quickPredict(nox.fit)
## plot result with method plot for quickPredict objects
plot(comp.pred)
## get prediction for the second term, equi using the binary version
equi.pred <- nox.fit %qp% "equi"
plot(equi.pred)
## get prediction for the interaction term, comp:equi
comp.equi.pred <- nox.fit %qp% "comp:equi"
## use image method image
image(comp.equi.pred)
## use contour method
contour(comp.equi.pred,col=2,lwd=2,labcex=1.5)
contour(comp.equi.pred,what="sd",lty=3,labcex=1.2,add=TRUE)
## use persp method
persp(comp.equi.pred,theta=-10,phi=20)
}
}
\keyword{hplot}
\keyword{methods}
\keyword{models}
\keyword{smooth}

