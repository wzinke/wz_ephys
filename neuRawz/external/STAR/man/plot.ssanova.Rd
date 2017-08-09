\name{plot.ssanova}
\Rdversion{1.1}
\alias{plot.ssanova}
\alias{plot.ssanova0}
\title{
A Plot Method for ssanova and ssanvoa0 Objects Tailored to Their Use in STAR
}
\description{
Plot a ssanova or a ssanova0 object.
}
\usage{
\method{plot}{ssanova}(x, y, include, ask = FALSE, ncol = 2, nrow = 3, ...)
\method{plot}{ssanova0}(x, y, include, ask = FALSE, ncol = 2, nrow = 3, ...)
}
\arguments{
  \item{x}{ a \code{\link[gss]{ssanova}} or a
    \code{\link[gss]{ssanova0}} object.}
  \item{y}{ not used, only included for compatibility with generic method.}
  \item{include}{ a character string with the model terms one wants to
    plot. If missing all terms are plotted.}
  \item{ask}{ a logical. If TRUE terms are plotted (on a common y scale)
    one after the other and the user is invited to hit the \code{enter}
    key to generate the next plot. If \code{FALSE} (default) all terms
    are drawn on a suitable number of X11 devices. The number of
    terms on each device is controlled by arguments ncol and nrow.}
  \item{ncol}{ the number of columns of the display matrix used on each device
    when \code{ask} is set to \code{FALSE}.}
  \item{nrow}{ the number of rows of the display matrix used on each device
    when \code{ask} is set to \code{FALSE}.}
  \item{\dots}{ not used only there for method definition compatibility.}
}
\value{
Nothing returned. The method is used for its side effect, plots are generated.
}
\author{
Christophe Pouzat \email{christophe.pouzat@gmail.com}
}
\note{
The designed is inspired by the \code{plot} method for \code{gam}
objects in package \code{mgcv}.
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
plot(m1.fit,nr=3,nc=1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
