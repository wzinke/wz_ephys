\name{quickPredict}
\Rdversion{1.1}
\alias{quickPredict}
\alias{\%qp\%}

\title{
A Simple Interface to predict method for ssanova and ssanova0 objects
}
\description{
Designed to quickly compute the effect of a \emph{single} model term.
This term can correspond to a single variable effect or to the
interaction of two variables.
}
\usage{
quickPredict(object, include = object$terms$labels[2], se.fit = TRUE,
             length.out, otherTermsFct = median)
object \%qp\% include
}

\arguments{
  \item{object}{ an object inheriting from ssanova and ssanova0
  (\code{\link[gss]{gssanova}} and \code{\link[gss]{gssanova0}} objects are therefore suitable).
}
  \item{include}{ a character string corresponding to a \emph{single}
  model term. See \code{\link[gss]{predict.ssanova}} and \code{\link[gss]{predict.ssanova}}. 
}
  \item{se.fit}{ logical flag indicating if standard errors are
  required. See \code{\link[gss]{predict.ssanova}} and \code{\link[gss]{predict.ssanova}}.
}
  \item{length.out}{ a positive integer, the number of points at which
  the prediction should be performed. These points are uniformly spread
  on the definition domain of the variable(s) implicitely specified by
  argument \code{include}. If \code{missing} a default of 501 for terms
  involving a single variable and of 101 for interaction terms involving
  two variables is provided.
}
  \item{otherTermsFct}{ a function applied to the other variables
  required for model specification.
}
}

\details{
\code{\%qp\%} is the binary version of \code{quickPredict}.
}

\value{
A \code{quickPredict} object. This object is a \code{\link{list}} with the following
components:
\item{xx }{a numeric vector with the values of the variable specified by
  the model term selected by argument \code{include}. When an interaction term was
  selected the values of the first variable are stored here.}
\item{yy }{a numeric vector with the values of the \emph{second} variable specified by
  the \emph{interaction} term selected by argument \code{include}. When
  selected term is not an interaction term, this component is
  \code{NULL}.}
\item{include }{ the value of the argument with this name.}
\item{call }{ the matched call.}
\item{est.mean }{ a numeric vector or matrix, for intercation terms,
  containing the estimated mean of the term.}
\item{est.sd }{ a numeric vector or matrix, for intercation terms,
  containing the estimated SD of the term. Is \code{NULL} is argument
  \code{se.fit} was \code{FALSE}.}
}
\author{
Christophe Pouzat \email{christophe.pouzat@gmail.com}
}
\seealso{
  \code{\link[gss]{predict.ssanova}},
  \code{\link[gss]{predict.ssanova}},
  \code{\link{plot.quickPredict}},
  \code{\link{image.quickPredict}},
  \code{\link{contour.quickPredict}},
  \code{\link{persp.quickPredict}},
  \code{\link{plot.ssanova}}
}
\examples{
## Follow up of ssanova example of gss
data(nox)
nox.fit <- ssanova(log10(nox)~comp*equi,data=nox)
## get prediction for the first term, comp
comp.pred <- quickPredict(nox.fit)
## plot result with method plot for quickPredict objects
plot(comp.pred)
## get prediction for the second term, equi using the binary version
equi.pred <- nox.fit \%qp\% "equi"
plot(equi.pred)
## get prediction for the interaction term, comp:equi
comp.equi.pred <- nox.fit \%qp\% "comp:equi"
## use image method image
image(comp.equi.pred)
## use contour method
contour(comp.equi.pred,col=2,lwd=2,labcex=1.5)
contour(comp.equi.pred,what="sd",lty=3,labcex=1.2,add=TRUE)
## use persp method
persp(comp.equi.pred,theta=-10,phi=20)
}
\keyword{models}
\keyword{smooth}
\keyword{regression}
