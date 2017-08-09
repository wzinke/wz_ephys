\name{predictLogProb}
\Rdversion{1.1}
\alias{predictLogProb}
\title{
Compute the Log Probability of a "New" Data Set
Using a Fitted Model Prediction
}
\description{
This function is designed to select models by cross-validation.
If two models A and B managed to pass the Ogata's tests, one should split the
data set in two parts, part1 and part2 fit each model on each part to get four
fitted model objects: A1, A2, B1, B2. The chosen model should then be the one
giving the largest of: \code{predictLogProb(A1,part2) + predictLogProb(A2,part1)} and
\code{predictLogProb(B1,part2) + predictLogProb(B2,part1)}.
}
\usage{
predictLogProb(object, newdata)
}
\arguments{
  \item{object}{ an object inheriting from \code{\link[gss]{ssanova}}
    and \code{\link[gss]{ssanova0}} (\code{\link[gss]{gssanova}} and
    \code{\link[gss]{gssanova0}} objects are therefore suitable).} 
  \item{newdata}{ a data frame containing the required variables. This
    data frame \emph{must be different from the one used to obtain} \code{object}.}
}
\details{
 If eta[i] is the prediction of the fitted model for element i of
 \code{newdata} the log probability is given by :
 
 \code{event[i] * eta[i] - log(1 + exp(eta[i]))}

 Where \code{event[i]} is 0 or 1 depending on the absence or presence of
 a spike at the considered time. 
 \emph{A binomial regression is assumed here}.
}
\value{
  A numeric, the sum over the index \code{i} above, the log probability
  of the data contained in \code{newdata} 
  assuming that the model contained in \code{object} is correct.
}
\author{
Christophe Pouzat \email{christophe.pouzat@gmail.com}
}

\seealso{
  \code{\link[STAR]{mkGLMdf}},
  \code{\link[STAR]{quickPredict}}
  \code{\link[gss]{gssanova}}
}
\examples{
\dontrun{
data(e060824spont)
summary(e060824spont[["neuron 1"]])
reportHTML(e060824spont[["neuron 1"]],filename="e060824spont_1",otherST=e060824spont[c(2)],maxiter=100)
acf(diff(e060824spont[["neuron 1"]]),type="partial")
DFA <- subset(mkGLMdf(e060824spont,0.004,0,59),neuron==1)
DFA <- within(DFA,i1 <- isi(DFA,lag=1))
DFA <- DFA[complete.cases(DFA),]
m2u1 <- mkM2U(DFA,"lN.1",0,29,seed=20061001)
m2ui <- mkM2U(DFA,"i1",0,29,maxiter=200,seed=20061001)
DFA <- within(DFA,e1t <- m2u1(lN.1))
DFA <- within(DFA,i1t <- m2ui(i1))
with(DFA,plot(ecdf(e1t),pch="."))
with(DFA,plot(ecdf(i1t),pch="."))
DFAts <- as.ts(apply(DFA[,c("e1t","i1t")],2,qnorm))
plot(filter(DFAts,rep(1/125,125)))
system.time(GF1 <- gssanova(event ~ e1t+i1t, data=subset(DFA,time<=29),family="binomial",seed=20061001))
tt.1 <- GF1 \%tt\% subset(DFA,time>29)
tt.1.summary <- summary(tt.1)
tt.1.summary
plot(tt.1.summary,which=c(1,2,4,6))
renewalTestPlot(tt.1$ppspFct())
plot(GF1,nc=1,nr=2)
system.time(GF2 <- gssanova(event ~ e1t+i1t, data=subset(DFA,time>29),family="binomial",seed=20061001))
tt.2 <- GF2 \%tt\% subset(DFA,time<=29)
tt.2.summary <- summary(tt.2)
tt.2.summary
plot(tt.2.summary,which=c(1,2,4,6))
renewalTestPlot(tt.2$ppspFct())
plot(GF2,nc=1,nr=2)
system.time(GF3 <- gssanova(event ~ e1t*i1t, data=subset(DFA,time<=29),family="binomial",seed=20061001))
tt.3 <- GF3 \%tt\% subset(DFA,time>29)
(tt.3.summary <- summary(tt.3))
plot(tt.3.summary,which=c(1,2,4,6))
renewalTestPlot(tt.3$ppspFct())
plot(GF3,nc=1,nr=3)
system.time(GF4 <- gssanova(event ~ e1t*i1t, data=subset(DFA,time>29),family="binomial",seed=20061001))
tt.4 <- GF4 \%tt\% subset(DFA,time<=29)
(tt.4.summary <- summary(tt.4))
plot(tt.4.summary,which=c(1,2,4,6))
renewalTestPlot(tt.4$ppspFct())
plot(GF4,nc=1,nr=3)
## Get the log probability of the data with the additive model
predictLogProb(GF1,newdata=subset(DFA,time>29))+predictLogProb(GF2,newdata=subset(DFA,time<=29))
## Get the log probability of the data with the non-additive model
predictLogProb(GF3,newdata=subset(DFA,time>29))+predictLogProb(GF4,newdata=subset(DFA,time<=29))
## The non additive model is the "best" so refit it to the whole data set
system.time(GF5 <- gssanova(event ~ e1t*i1t, data=DFA,family="binomial",seed=20061001))
plot(GF5,nr=3,nc=1)
}
}
\keyword{models}
\keyword{regression}
