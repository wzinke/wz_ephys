\name{gssObj}
\alias{gssObj}
\alias{gssObj.gsspsth}
\alias{gssObj.gsspsth0}
\alias{gssObj.gsslockedTrain}
\alias{gssObj.gsslockedTrain0}
\title{Generic Function and Methods for Extracting a gss object}
\description{
  Some functions of \code{STAR}, like \code{gsspsth},  \code{gsspsth0} and
  \code{gsslockedTrain}, \code{gsslockedTrain0} perform gss fits internally and keep as a list
  component or within the environment of a returned function the result
  of this fit. Method \code{gssObj} extracts this \code{gss} object.
}
\usage{
gssObj(object, ...)
\method{gssObj}{gsspsth}(object, ...)
\method{gssObj}{gsspsth0}(object, ...)
\method{gssObj}{gsslockedTrain}(object, ...)
\method{gssObj}{gsslockedTrain0}(object, ...)
}
\arguments{
  \item{object}{an object containing a \code{gssanova} or a
    \code{gssanova0} object. Currently the
    result of a call to \code{\link{gsspsth}}, \code{\link{gsspsth0}} or to
    \code{\link{gsslockedTrain}}, \code{\link{gsslockedTrain0}}.
  }
  \item{\dots}{not used for now}
}
\value{
  A \code{\link[gss]{gssanova}} or a \code{\link[gss]{gssanova0}} object. 
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com}}
\seealso{
  \code{\link[gss]{gssanova}},
  \code{\link[gss]{gssanova0}},
  \code{\link{gsspsth}},
  \code{\link{gsspsth0}},
  \code{\link{gsslockedTrain}},
  \code{\link{gsslockedTrain0}}
}
\examples{
##
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{smooth}
\keyword{regression}
