\name{mkAR}
\Rdversion{1.1}
\alias{mkAR}
\title{
  Generate a Data Frame With Variables Suitable For an AR Like Model
}
\description{
  The variables added to the data frame corresponding to the first
  argument of the function are the former inter spike intervals. These variables
  are moreover transformed with \code{mkM2U} so that they have an approximately uniform
  distribution on their definition domain.
}
\usage{
mkAR(df, low, high, max.order, selfName = "lN.1",...)
}
\details{
  When \code{max.order} > 1 the previous inter spike intervals are all
  transformed using the "map to uniform" function estimated from the
  inter spike intervals at lag 1.  
}
\arguments{
  \item{df}{ a data frame. This data frame should contain a variable
    \code{time} like data frames returned by
    \code{\link[STAR]{mkGLMdf}}.
  }
  \item{low}{ a numeric, the smallest value of variable \code{time} from
  which the transformation is looked for. If missing defaults to the
  smallest time.}
  \item{high}{ a numeric, the largest value of variable \code{time} up to
  which the transformation is looked for. If missing defaults to the
  largest time.}
  \item{max.order}{ a postive integer, the maximal order of the AR
  model. How many previous inter spike intervals should be used in order
  to predict the duration of the next interval?}
  \item{selfName}{ a character string or an integer specifying the
  variable of \code{df} containing the elapsed time since the last spike
  of the considered neuron.}
  \item{\dots}{ additional arguments passed to \code{\link[STAR]{mkM2U}}}.}
\value{
  A data frame is returned. In addition to the variables of df the returned data frame contains
  a variable \code{est} with the transformed elapsed time since the last
  spike of the neuron and \code{i1t}, \code{i2t},\dots,\code{i max.order t}, the
  transformed previous inter spike intervals.
  The returned data frame has also four attributes:
  \item{fmla}{ a formula suitable for a first argument of, say, \code{\link[gss]{gssanova}}.}
  \item{m2uL}{ the function returned by mkM2U transforming the elasped time
    since the last spike of the neuron.}
  \item{m2uI}{ the function returned by mkM2U transforming the first former
    inter spike interval.}
  \item{call}{ the matched call.}
}
\author{
Christophe Pouzat \email{christophe.pouzat@gmail.com}
}

\seealso{
  \code{\link[STAR]{mkM2U}},
  \code{\link[gss]{gssanova}}
}
\examples{
\dontrun{
require(STAR)
data(e060824spont)
DFA <- subset(mkGLMdf(e060824spont,0.004,0,59),neuron==1)
DFA <- mkAR(DFA, 0, 29, 5, maxiter=200)
head(DFA)
tail(DFA)
ar.fit <- gssanova(attr(DFA,"fmla"), data=DFA,family="binomial",seed=20061001)
plot(ar.fit \%qp\% "est")
plot(ar.fit \%qp\% "i1t")
plot(ar.fit \%qp\% "i2t")
plot(ar.fit \%qp\% "i3t")
plot(ar.fit \%qp\% "i4t")
plot(ar.fit \%qp\% "i5t")
}
}
\keyword{models}
\keyword{regression}
\keyword{ts}
