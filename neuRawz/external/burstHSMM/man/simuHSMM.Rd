\name{simuHSMM} 
\Rdversion{0.1-1} 
\alias{simu.hsmm}
\title{Simulate spike-trains from Hidden Semi-Markov Model} 
\description{ 
To simulate spike trains from a hidden semi-Markov model
} 
\usage{ 
simu.hsmm(total.time, pars = NULL)
} 
\arguments{ 
\item{total.time}{a positive number indicating the time period in seconds.}
\item{pars}{parameters of the HSMM. See fit.hsmm for more details.}
}
\seealso{ 
\code{\link{fit.hsmm}}. 
} 
\value{A list containing
\item{spikes}{a numeric vector of the ISIs measured in ms.}
\item{states}{a binary vector burst/non-burst states of the ISIs.}
\item{transitions}{a numeric vector showing the state transition times.}
}
\references{Tokdar, S. T., Xi, P., Kelly, R. C. and Kass, R. E. (2009). Detection of bursts in extracellular spike trains using hidden semi-Markov point process models. \emph{Journal of Computational Neuroscience.}  DOI: 10.1007/s10827-009-0182-2} 
\examples{ 
## Generate a spike over 10s
## burst mean / SD in sec: (Burst) 0.2 / 0.04; (Non-burst) 0.7 / 0.07
## ISI mean / SD in sec: (Burst) 0.02 / 0.01; (Non-burst) 0.10 / 0.10
spk <- simu.hsmm(10, pars = c(25, 0.008, 100, 0.007, 4, 0.005, 1, .1))

## plot
time.stamps <- cumsum(spk$spikes)
plot(stepfun(time.stamps, c(spk$states, spk$states[length(spk$states)])), do.points = FALSE)
} 
