\name{cockroachAlData}
\alias{CAL1S}
\alias{CAL1V}
\alias{CAL2S}
\alias{CAL2C}
\alias{e060517spont}
\alias{e060517ionon}
\alias{e060817spont}
\alias{e060817terpi}
\alias{e060817citron}
\alias{e060817mix}
\alias{e060824spont}
\alias{e060824citral}
\alias{e070528spont}
\alias{e070528citronellal}
\docType{data}
\title{Spike Trains of several Cockroach Antennal Lobe Neurons
  Recorded from Six Animals}
\description{
  Four (\code{CAL1S} and \code{CAL1V}), three (\code{CAL2S} and
  \code{CAL2C}), three (\code{e060517spont} and \code{e060517ionon}),
  three (\code{e060817spont}, \code{e060817terpi}, \code{e060817citron}
  and \code{e060817mix}), two (\code{e060824spont} and
  \code{e060824citral}) and four (\code{e070528spont} and
  \code{e070528citronellal}) Cockroach (\emph{Periplaneta americana}) antennal lobe neurons
  (putative projection neurons) were recorded simultaneously
  and extracellularly during spontaneous activity and odors (vanilin,
  citral, citronellal, terpineol, beta-ionon)
  responses from six different animals. The data sets contain the sorted
  spike trains of the neurons. 
}
\usage{
data(CAL1S)
data(CAL1V)
data(CAL2S)
data(CAL2C)
data(e060517spont)
data(e060517ionon)
data(e060817spont)
data(e060817terpi)
data(e060817citron)
data(e060817mix)
data(e060824spont)
data(e060824citral)
data(e070528spont)
data(e070528citronellal)
}
\format{
  \code{CAL1S} is a named list with 4 components
  (\code{"neuron 1"}, \code{"neuron 2"}, \code{"neuron 3"},
  \code{"neuron 4"}). Each
  component contains the spike train (ie, action potentials occurrence
  times) of one neuron recorded during 30 s of spontaneous
  activity. \emph{Times are expressed in seconds}.

  \code{CAL1V} is a named list with 4 components
  (\code{"neuron 1"}, \code{"neuron 2"}, \code{"neuron 3"},
  \code{"neuron 4"}).
  Each component is a named list
  with 20 components: \code{"stim. 1"}, ..., \code{"stim. 20"}. Each
  sub-list contains the spike train of one neuron during 1 stimulation
  (odor puff) with \emph{vanillin}
  (\url{http://en.wikipedia.org/wiki/Vanillin}). Each acquisition was 10
  s long. The command 
  to the odor delivery  valve was on between sec 4.49 and sec 4.99. 

  \code{CAL2S} is a named list with 3 components
  (\code{"neuron 1"}, \code{"neuron 2"}, \code{"neuron 3"}). Each
  component contains the spike train (ie, action potentials occurrence
  times) of one neuron recorded during 1 mn of spontaneous
  activity. \emph{Times are expressed in seconds}.

  \code{CAL2C} is a named list with 3 components
  (\code{"neuron 1"}, \code{"neuron 2"}, \code{"neuron 3"}).
  Each component is a named list
  with 20 components: \code{"stim. 1"}, ..., \code{"stim. 20"}. Each
  sub-list contains the spike train of one neuron during 1 stimulation
  (odor puff) with \emph{citral}
  (\url{http://en.wikipedia.org/wiki/Citral}). Each acquisition was 14 s long. The command
  to the odor delivery  valve was on between sec 5.87 and sec 6.37.
  
  \code{e060517spont} is a named list of  with 3 components
  (\code{"neuron 1"}, \code{"neuron 2"}, \code{"neuron 3"}). Each
  component is a \code{spikeTrain} object (ie, action potentials occurrence
  times) of one neuron recorded during 61 s of spontaneous
  activity. \emph{Times are expressed in seconds}.

  \code{e060517ionon} is a named list with 3 components
  (\code{"neuron 1"}, \code{"neuron 2"}, \code{"neuron 3"}).
  Each component is a \code{repeatedTrain} object
  with 19 \code{spikeTrain} objects: \code{"stim. 1"}, ..., \code{"stim. 19"}. Each
  \code{spikeTrain} contains the spike train of one neuron during 1 stimulation
  (odor puff) with \emph{beta-ionon}
  (\url{http://commons.wikimedia.org/wiki/Image:Beta-Ionon.svg}). Each acquisition was 15 s long. The command
  to the odor delivery  valve was on between sec 6.07 and sec 6.57.

  \code{e060817spont} is a named list of  with 3 components
  (\code{"neuron 1"}, \code{"neuron 2"}, \code{"neuron 3"}). Each
  component is a \code{spikeTrain} object (ie, action potentials occurrence
  times) of one neuron recorded during 60 s of spontaneous
  activity. \emph{Times are expressed in seconds}.

  \code{e060817terpi} is a named list with 3 components
  (\code{"neuron 1"}, \code{"neuron 2"}, \code{"neuron 3"}).
  Each component is a \code{repeatedTrain} object
  with 20 \code{spikeTrain} objects: \code{"stim. 1"}, ..., \code{"stim. 20"}. Each
  \code{spikeTrain} contains the spike train of one neuron during 1 stimulation
  (odor puff) with \emph{terpineol}
  (\url{http://en.wikipedia.org/wiki/Terpineol}). Each acquisition was 15 s long. The command
  to the odor delivery valve was on between sec 6.03 and sec 6.53.

  \code{e060817citron} is a named list with 3 components
  (\code{"neuron 1"}, \code{"neuron 2"}, \code{"neuron 3"}).
  Each component is a \code{repeatedTrain} object
  with 20 \code{spikeTrain} objects: \code{"stim. 1"}, ..., \code{"stim. 20"}. Each
  \code{spikeTrain} contains the spike train of one neuron during 1 stimulation
  (odor puff) with \emph{citronellal}
  (\url{http://en.wikipedia.org/wiki/Citronellal}). Each acquisition was 15 s long. The command
  to the odor delivery valve was on between sec 5.99 and sec 6.49.

  \code{e060817mix} is a named list with 3 components
  (\code{"neuron 1"}, \code{"neuron 2"}, \code{"neuron 3"}).
  Each component is a \code{repeatedTrain} object
  with 20 \code{spikeTrain} objects: \code{"stim. 1"}, ..., \code{"stim. 20"}. Each
  \code{spikeTrain} contains the spike train of one neuron during 1 stimulation
  (odor puff) with a mixture of \emph{terpinaol} and \emph{citronellal}
  (the sum of the two previous stim.). Each acquisition was 15 s long. The command
  to the odor delivery valve was on between sec 6.01 and sec 6.51.

  \code{e060824spont} is a named list of  with 2 components
  (\code{"neuron 1"}, \code{"neuron 2"}). Each
  component is a \code{spikeTrain} object (ie, action potentials occurrence
  times) of one neuron recorded during 59 s of spontaneous
  activity. \emph{Times are expressed in seconds}.

  \code{e060824citral} is a named list with 2 components
  (\code{"neuron 1"}, \code{"neuron 2"}).
  Each component is a named list
  with 20 components: \code{"stim. 1"}, ..., \code{"stim. 20"}. Each
  sub-list contains the spike train of one neuron during 1 stimulation
  (odor puff) with \emph{citral}
  (\url{http://en.wikipedia.org/wiki/Citral}). Each acquisition was 15 s long. The command
  to the odor delivery  valve was on between sec 6.01 and sec 6.51.
  
  \code{e070528spont} is a named list of  with 4 components
  (\code{"neuron 1"}, \code{"neuron 2"}, \code{"neuron 3"},
  \code{"neuron 4"}). Each
  component is a \code{spikeTrain} object (ie, action potentials occurrence
  times) of one neuron recorded during 60 s of spontaneous
  activity. \emph{Times are expressed in seconds}.

  \code{e070528citronellal} is a named list with 4 components
  (\code{"neuron 1"}, \code{"neuron 2"}, \code{"neuron 3"},
  \code{"neuron 4"}).
  Each component is a \code{repeatedTrain} object
  with 15 \code{spikeTrain} objects: \code{"stim. 1"}, ..., \code{"stim. 15"}. Each
  \code{spikeTrain} contains the spike train of one neuron during 1 stimulation
  (odor puff) with \emph{citronellal}
  (\url{http://en.wikipedia.org/wiki/Citronellal}). Each acquisition was 13 s long. The command
  to the odor delivery  valve was on between sec 6.14 and sec 6.64.
}
\details{
  Every \code{repeatedTrain} object of these data sets has an \code{attribute} named
  \code{stimTimeCourse} containing the openng and closing times of the
  odor delivery valve.
  
  The data were recorded from neighboring sites on a \emph{NeuroNexus}
  (\url{http://neuronexustech.com/}) silicon probe. Sorting was done
  with \code{SpikeOMatic} with superposition resolution which can AND
  DOES lead to artifcats on cross-correlograms. 
}
\source{
  Recording and spike sorting performed by Antoine Chaffiol
  \email{antoine.chaffiol@univ-paris5.fr} at the Cerebral Physiology
  Lab, CNRS UMR 8118:
  \url{http://www.biomedicale.univ-paris5.fr/physcerv/physiologie_cerebrale.htm}. 
}
\references{
  \url{http://www.biomedicale.univ-paris5.fr/physcerv/C_Pouzat/Doc/ChaffiolEtAl_FENS2006.pdf}
}
\examples{
## load CAL1S data
data(CAL1S)
## convert the data into spikeTrain objects
CAL1S <- lapply(CAL1S,as.spikeTrain)
## look at the train of the 1sd neuron
CAL1S[["neuron 1"]]
## fit the 6 different renewal models to the 1st neuron spike train
compModels(CAL1S[["neuron 1"]])
## look at the ISI distribution with the fitted invgauss dist for
## this 1st neuron
isiHistFit(CAL1S[["neuron 1"]],model="invgauss")

## load CAL1V data
data(CAL1V)
## convert them to repeatedTrain objects
CAL1V <- lapply(CAL1V, as.repeatedTrain)
## look at the raster of the 1st neuron
CAL1V[["neuron 1"]]

## load e070528spont data
data(e070528spont)
## look at the spike train of the 1st neuron
e070528spont[["neuron 1"]]

## load e070528citronellal data
data(e070528citronellal)
## Get the stimulus time course
attr(e070528citronellal[["neuron 1"]],"stimTimeCourse")
## look at the raster of the 1st neuron
plot(e070528citronellal[["neuron 1"]],stim=c(6.14,6.64))

\dontrun{
## A "detailed" analysis of e060817 were 2 odors as well as there mixtures
## were used.
## Load the terpineol, citronellal and mixture response data
data(e060817terpi)
data(e060817citron)
data(e060817mix)
## get smooth psths with gsspsth0
e060817terpiN1PSTH <- gsspsth0(e060817terpi[["neuron 1"]])
e060817terpiN2PSTH <- gsspsth0(e060817terpi[["neuron 2"]])
e060817terpiN3PSTH <- gsspsth0(e060817terpi[["neuron 3"]])
e060817citronN1PSTH <- gsspsth0(e060817citron[["neuron 1"]])
e060817citronN2PSTH <- gsspsth0(e060817citron[["neuron 2"]])
e060817citronN3PSTH <- gsspsth0(e060817citron[["neuron 3"]])
e060817mixN1PSTH <- gsspsth0(e060817mix[["neuron 1"]])
e060817mixN2PSTH <- gsspsth0(e060817mix[["neuron 2"]])
e060817mixN3PSTH <- gsspsth0(e060817mix[["neuron 3"]])
## look at them
## Neuron 1
plot(e060817terpiN1PSTH,stimTimeCourse=attr(e060817terpi[["neuron 1"]],"stimTimeCourse"),colCI=2)
plot(e060817citronN1PSTH,stimTimeCourse=attr(e060817citron[["neuron 1"]],"stimTimeCourse"),colCI=2)
plot(e060817mixN1PSTH,stimTimeCourse=attr(e060817mix[["neuron 1"]],"stimTimeCourse"),colCI=2)
## Neuron 2
plot(e060817terpiN2PSTH,stimTimeCourse=attr(e060817terpi[["neuron 2"]],"stimTimeCourse"),colCI=2)
plot(e060817citronN2PSTH,stimTimeCourse=attr(e060817citron[["neuron 2"]],"stimTimeCourse"),colCI=2)
plot(e060817mixN2PSTH,stimTimeCourse=attr(e060817mix[["neuron 2"]],"stimTimeCourse"),colCI=2)
## Neuron 3
plot(e060817terpiN3PSTH,stimTimeCourse=attr(e060817terpi[["neuron 3"]],"stimTimeCourse"),colCI=2)
plot(e060817citronN3PSTH,stimTimeCourse=attr(e060817citron[["neuron 3"]],"stimTimeCourse"),colCI=2)
plot(e060817mixN3PSTH,stimTimeCourse=attr(e060817mix[["neuron 3"]],"stimTimeCourse"),colCI=2)

## Make now fancier plots with superposed psths ####
## Take into account the fact that the stimuli onsets are not identical

## Neuron 1
plot(e060817mixN1PSTH$mids-0.02,e060817mixN1PSTH$ciUp,type="n",ylim=c(0,max(e060817mixN1PSTH$ciUp)),xlim=c(5,14),xlab="Time (s)",ylab="Firing rate (Hz)",main="Neuron 1 e060817")
rect(5.99,0,6.49,max(e060817mixN1PSTH$ciUp),col="grey80",border=NA)
abline(h=0)
polygon(c(e060817mixN1PSTH$mids-0.02,rev(e060817mixN1PSTH$mids-0.02)),c(e060817mixN1PSTH$ciLow,rev(e060817mixN1PSTH$ciUp)),col=rgb(1,0,1,0.5),border=NA)
polygon(c(e060817citronN1PSTH$mids,rev(e060817citronN1PSTH$mids)),c(e060817citronN1PSTH$ciLow,rev(e060817citronN1PSTH$ciUp)),col=rgb(1,0,0,0.5),border=NA)
polygon(c(e060817terpiN1PSTH$mids-0.04,rev(e060817terpiN1PSTH$mids-0.04)),c(e060817terpiN1PSTH$ciLow,rev(e060817terpiN1PSTH$ciUp)),col=rgb(0,0,1,0.5),border=NA)
lines(e060817terpiN1PSTH$mids-0.04,e060817terpiN1PSTH$freq,col=rgb(0,0,1),lwd=2)
lines(e060817citronN1PSTH$mids,e060817citronN1PSTH$freq,col=rgb(1,0,0),lwd=2)
lines(e060817mixN1PSTH$mids-0.02,e060817mixN1PSTH$freq,col=rgb(0,0,0),lwd=2)
legend(8,0.9*max(e060817mixN1PSTH$ciUp),c("Terpineol","Citronellal","Mixture"),col=c(4,2,1),lwd=2)

## Neuron 2
plot(e060817mixN2PSTH$mids-0.02,e060817mixN2PSTH$ciUp,type="n",ylim=c(0,max(e060817mixN2PSTH$ciUp)),xlim=c(5,14),xlab="Time (s)",ylab="Firing rate (Hz)",main="Neuron 2 e060817")
rect(5.99,0,6.49,max(e060817mixN2PSTH$ciUp),col="grey80",border=NA)
abline(h=0)
polygon(c(e060817mixN2PSTH$mids-0.02,rev(e060817mixN2PSTH$mids-0.02)),c(e060817mixN2PSTH$ciLow,rev(e060817mixN2PSTH$ciUp)),col=rgb(1,0,1,0.5),border=NA)
polygon(c(e060817citronN2PSTH$mids,rev(e060817citronN2PSTH$mids)),c(e060817citronN2PSTH$ciLow,rev(e060817citronN2PSTH$ciUp)),col=rgb(1,0,0,0.5),border=NA)
polygon(c(e060817terpiN2PSTH$mids-0.04,rev(e060817terpiN2PSTH$mids-0.04)),c(e060817terpiN2PSTH$ciLow,rev(e060817terpiN2PSTH$ciUp)),col=rgb(0,0,1,0.5),border=NA)
lines(e060817terpiN2PSTH$mids-0.04,e060817terpiN2PSTH$freq,col=rgb(0,0,1),lwd=2)
lines(e060817citronN2PSTH$mids,e060817citronN2PSTH$freq,col=rgb(1,0,0),lwd=2)
lines(e060817mixN2PSTH$mids-0.02,e060817mixN2PSTH$freq,col=rgb(0,0,0),lwd=2)
legend(8,0.9*max(e060817mixN2PSTH$ciUp),c("Terpineol","Citronellal","Mixture"),col=c(4,2,1),lwd=2)

## Neuron 3
plot(e060817mixN3PSTH$mids-0.02,e060817mixN3PSTH$ciUp,type="n",ylim=c(0,max(e060817mixN3PSTH$ciUp)),xlim=c(5,14),xlab="Time (s)",ylab="Firing rate (Hz)",main="Neuron 3 e060817")
rect(5.99,0,6.49,max(e060817mixN3PSTH$ciUp),col="grey80",border=NA)
abline(h=0)
polygon(c(e060817mixN3PSTH$mids-0.02,rev(e060817mixN3PSTH$mids-0.02)),c(e060817mixN3PSTH$ciLow,rev(e060817mixN3PSTH$ciUp)),col=rgb(1,0,1,0.5),border=NA)
polygon(c(e060817citronN3PSTH$mids,rev(e060817citronN3PSTH$mids)),c(e060817citronN3PSTH$ciLow,rev(e060817citronN3PSTH$ciUp)),col=rgb(1,0,0,0.5),border=NA)
polygon(c(e060817terpiN3PSTH$mids-0.04,rev(e060817terpiN3PSTH$mids-0.04)),c(e060817terpiN3PSTH$ciLow,rev(e060817terpiN3PSTH$ciUp)),col=rgb(0,0,1,0.5),border=NA)
lines(e060817terpiN3PSTH$mids-0.04,e060817terpiN3PSTH$freq,col=rgb(0,0,1),lwd=2)
lines(e060817citronN3PSTH$mids,e060817citronN3PSTH$freq,col=rgb(1,0,0),lwd=2)
lines(e060817mixN3PSTH$mids-0.02,e060817mixN3PSTH$freq,col=rgb(0,0,0),lwd=2)
legend(8,0.9*max(e060817mixN3PSTH$ciUp),c("Terpineol","Citronellal","Mixture"),col=c(4,2,1),lwd=2)
}

}
\keyword{datasets}
