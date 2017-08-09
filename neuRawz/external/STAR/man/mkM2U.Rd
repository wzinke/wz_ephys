\name{mkM2U}
\Rdversion{1.1}
\alias{mkM2U}
\title{
Makes a Smooth Function Mapping a Data Frame Variable Onto a Variable
Uniform on Its Definition Domain 
}
\description{
The smooth transformation function is a smooth version of the
\code{\link{ecdf}}. A smooth density estimate as well as the inverse
transformation (the quantile function) are also returned as attributes. 
}
\usage{
mkM2U(df, vN, low, high, delta, alpha=2, ...)
}

\arguments{
  \item{df}{ a data frame. This data frame should contain a variable
    \code{time} like data frames returned by \code{\link[STAR]{mkGLMdf}}.}
  \item{vN}{ a character string corresponding to the name of one of the
  variables of \code{df} or an integer, its index. Variable \code{vN} is
  the one for which the mapping to uniform is looked for.}
  \item{low}{ a numeric, the smallest value of variable \code{time} from
  which the transformation is looked for. If missing defaults to the
  smallest time.}
  \item{high}{ a numeric, the largest value of variable \code{time} up to
  which the transformation is looked for. If missing defaults to the
  largest time.}
  \item{delta}{ a numeric, the bin width used to build the variable
  values histogram. This histogram is subsequently smoothed. Default
  provided if missing.}
  \item{alpha}{ see \code{\link[gss]{ssden}}.}
  \item{\dots}{ additional arguments passed to \code{\link[gss]{ssden}} called
  internally by the function}
}
\details{
  The smooth mapping to uniform function returned by \code{mkM2U} is
  obtained by first selecting a subset of the variable values for which
  the variable \code{time} of \code{df} is between \code{low} and
  \code{high}. The values are then binned between the \code{min} and the
  \code{max} of the (complete) variable values with a bin width
  \code{delta}. Function \code{\link[gss]{ssden}} is then called on the
  histogram and the result is stored in object \code{ii.fit} (This
  object is stored in the \code{closure} of the returned function). The returned
  function is the result of a call of \code{pssden} on
  ii.fit and the argument.
  
  A function inverting the "mapping to uniform function", that is, a
  quantile function, is also returned as \link{attributes}
  \code{qFct}. This 
  inverse function is obtained by numerical inversion, calling
  \code{\link{uniroot}} internally. Additional arguments can be passed
  to \code{\link{uniroot}} via the \dots argument of the function.

  A function returning the smooth density estimate is returned as
  \link{attributes} \code{dFct}.
}
\value{
  A function returning the probability for \code{vN} random variable to
  have a value smaller or equal to its first argument. The returned
  function calls internally \code{\link{integrate}}. Additional
  arguments can be passed to the latter via the \dots argument of the
  returned function.

  As explained in the \code{details} section, the returned function has
  the smooth density function, \code{dFct}, as well as the inverse
  function, \code{qFct}, as attributes. Attribute \code{call} contains
  the matched call and \code{range} contains the full range of the
  mapped variable.
}
\author{
Christophe Pouzat \email{christophe.pouzat@gmail.com}
}
\note{
  Since the density returned by \code{\link[gss]{dssden}} can sometime
  integrate to a value slightly different from 1 on its definition
  domain, the actual integral is evaluated with \code{\link{integrate}}
  and the returned density is renormalised. A look-up table of 101
  regularly spaced quantiles and the corresponding probabilities is also
  created and stored in the returned function closure. This look-up
  table is used to speed up the computations performed by the returned
  function which uses \code{\link{integrate}} and not
  \code{pssden}. It is also used to speed up the
  computations of the inverse function (returned as attribute
  \code{qFct}) which uses \code{\link{uniroot}} and not
  \code{qssden}.
}

\seealso{
  \code{\link[gss]{ssden}},
  \code{\link[gss]{dssden}},
  \code{\link{integrate}},
  \code{\link{uniroot}},
  \code{\link[STAR]{mkGLMdf}}
}
\examples{
\dontrun{
require(STAR)
data(e060824spont)
DFA <- subset(mkGLMdf(e060824spont,0.004,0,59),neuron==1)
DFA <- within(DFA,i1 <- isi(DFA,lag=1))
DFA <- DFA[complete.cases(DFA),]
m2u1 <- mkM2U(DFA,"lN.1",0,29)
m2ui <- mkM2U(DFA,"i1",0,29,maxiter=200)
DFA <- within(DFA,e1t <- m2u1(lN.1))
DFA <- within(DFA,i1t <- m2ui(i1))
with(DFA,plot(ecdf(e1t[time>29]),pch="."))
abline(a=0,b=1,col=2,lty=2)
with(DFA,plot(ecdf(i1t[time>29]),pch="."))
abline(a=0,b=1,col=2,lty=2)}
}
\keyword{models}
\keyword{regression}

