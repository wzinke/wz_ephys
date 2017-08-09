\name{mkDummy}
\alias{mkDummy}
\title{Generates a Data Frame of Dummy Variables for Use in gam}
\description{
  Using argument \code{by} in \code{\link[mgcv]{s}} or
  \code{\link[mgcv]{te}} of \code{\link[mgcv]{gam}} requires dummy
  variables to be set up. This is the job of this function.
}
\usage{
mkDummy(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a \code{factor}.}
}
\value{
  A \code{\link{data.frame}} with as many variables as there are
  \code{levels} in \code{x} and as many rows as elements in \code{x}.
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com}}
\seealso{
  \code{\link{mkGLMdf}},
  \code{\link[mgcv]{gam}},
  \code{\link[mgcv]{s}},
  \code{\link[mgcv]{te}}
}
\examples{
## coming soon
}
\keyword{models}
