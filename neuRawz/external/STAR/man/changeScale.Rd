\name{changeScale}
\Rdversion{1.1}
\alias{changeScale}
\title{
Change the Scales of a quickPredict Object for an Interaction Term
}
\description{
  Designed to transform results of \code{\link[STAR]{quickPredict}} obtained on interaction
  terms from the transformed scale (on which the variables are
  approximately uniformly distributed) onto the "native", linear scale. 
}
\usage{
changeScale(obj, xFct, yFct)
}
\arguments{
  \item{obj}{ a \code{\link[STAR]{quickPredict}} object.}
  \item{xFct}{ a function to be applied on the \code{xx} element of
    \code{obj}. This function should be the \code{qFct} attribute of the
    function, returned by \code{\link[STAR]{mkM2U}}, used to transform
    the variable from the "native" to the "uniform" scale.}
  \item{yFct}{ a function to be applied on the \code{yy} element of
    \code{obj}. This function should be the \code{qFct} attribute of the
    function, returned by \code{\link[STAR]{mkM2U}}, used to transform
    the variable from the "native" to the "uniform" scale.}
}
\value{
  A \code{\link[STAR]{quickPredict}} object.
}
\author{
Christophe Pouzat \email{christophe.pouzat@gmail.com}
}

\seealso{
  \code{\link[STAR]{quickPredict}},
  \code{\link[STAR]{plot.quickPredict}}
}
\examples{
\dontrun{
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
abline(a=0,b=1,col=2,lty=2)
m1.fit <- gssanova(event~e1t*i1t, data=subset(DFA,time>29), family="binomial", seed=20061001)
inter.pred <- m1.fit \%qp\% "e1t:i1t"
contour(inter.pred,what="mean",nlevels=10,col=2,lwd=2)
contour(inter.pred,what="sd",nlevels=5,col=1,lwd=1,lty=2,add=TRUE)
inter.predN <- changeScale(inter.pred,attr(m2u1,"qFct"),attr(m2ui,"qFct"))
contour(inter.predN,what="mean",nlevels=5,col=2,lwd=1)
contour(inter.predN,what="sd",nlevels=3,col=1,lwd=1,lty=2,add=TRUE)
}
}
\keyword{models}
\keyword{regression}
