\name{purkinjeCellData}
\alias{sPK}
\alias{mPK}
\docType{data}
\title{Spike Trains of a Purkinje Cells (PC) Recorded in Control
  Conditions and With Bath Applied Bicuculline} 
\description{
    An object of class \code{"SpikeTrain"}. Spontaneous discharge of a
  single PC recorded during 300 s in normal saline conditions and during
  300 s in the presence of 25 \eqn{\mu}{micro}M bath applied bicuculline.
}
\usage{
data(sPK)
data(mPK)
}
\format{
  \code{sPK} is a named list with 2 components
  (\code{"ctl"}, \code{"bicu"}. Each
  component contains the spike train (ie, action potentials occurrence
  times) of one Purkinje cell recorded during 300 s of spontaneous
  activity in control (\code{"ctl"}) condition and with bath applied
  bicuculline (\code{"bicu"}). \emph{Times are expressed in seconds}.

  \code{mPK} is a named list with 8 components
  (\code{"neuron 1"}, \code{"neuron 2"}, ..., \code{"neuron 8"}. Each
  component is itself a list with  the spike train (ie, action potentials occurrence
  times) of one Purkinje cell recorded during 300 s of spontaneous
  activity in control (\code{"ctl"}) condition and with bath applied
  bicuculline (\code{"bicu"}). \emph{Times are expressed in seconds}.

}
\details{
  The recording contained in \code{sPK} was done in cell-attached
  mode. The one in \code{mPK} was done with a NeuroNexus silicon probe.
  
  Bicuculline is a GABAA receptor antagonist. It blocks all GABAA inhibition.
}
\source{
  Recording and spike sorting performed by Matthieu Delescluse
  at the Cerebral Physiology Lab, CNRS UMR 8118:
  \url{http://www.biomedicale.univ-paris5.fr/physcerv/physiologie_cerebrale.htm}. 
}
\examples{
\dontrun{
## load spontaneous data of 1 Purkinje cell
## recorded in cell attached mode from a cerebellar
## slice in control and bath applied bicuculline conditions
data(sPK)
## coerce data to spikeTrain objects
sPK <- lapply(sPK,as.spikeTrain)
## Get a summary of the ctl data
summary(sPK[["ctl"]])
## Look at the control train
## Don't show the rug plot for clarity
plot(sPK[["ctl"]],addRug=FALSE)
## Generate the renewal test plot taking into account
## the size of the data set (a lot of spikes!).
renewalTestPlot(sPK[["ctl"]],d=10,orderPlotPch=".",lag.max=250)
## Get a summary of the bicu data
summary(sPK[["bicu"]])
## Look at the control train
## Don't show the rug plot for clarity
plot(sPK[["bicu"]],addRug=FALSE)
## Generate the renewal test plot taking into account
## the size of the data set (a lot of spikes!).
renewalTestPlot(sPK[["bicu"]],d=10,orderPlotPch=".",lag.max=250);par(oldpar)
## This time the data are NOT stationary. This is seen clearly on a acf
## plot with very large lag.max
acf.spikeTrain(sPK[["bicu"]],lag.max=2000)

}
}
\keyword{datasets}

