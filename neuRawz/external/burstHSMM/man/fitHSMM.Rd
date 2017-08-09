\name{fitHSMM}
\Rdversion{0.1-1}
\alias{fit.hsmm}
\alias{print.hsmm}
\alias{plot.hsmm}
\alias{summary.hsmm} 
\title{Hidden Semi-Markov Model}
\description{ 
Fits a hidden semi-Markov model to a spike train via MCMC
} 
\usage{ 
fit.hsmm(y, scan = c(5, 25), dur.shape = c(1, 1), limit = NULL,
	states = NULL, pars = NULL, hpar = NULL, tune = NULL,
	nburn = 0, nsamp = 1e3, nskip = 1, nshot = 5)
\method{print}{hsmm}(obj)
\method{plot}{hsmm}(obj, subset = NULL, glen = 1e2, newfig = FALSE, hbreaks = NULL)
\method{summary}{hsmm}(obj, subset = NULL, newfig = FALSE, toplot = TRUE)
} 
\arguments{
\item{y}{a numeric vector containing the inter-Spike intervals.}
\item{scan}{a numeric vector specifying the minimum and maximum 
scan length.}
\item{dur.shape}{a numeric vector specifying the shape parameters for the 
  gamma distributions for modeling burst and non-burst durations.}
\item{limit}{a numeric value used as the maximum limit for a bursting ISI.}
\item{states}{a binary vector to initialize the hidden states for MCMC 
iterations.}
\item{pars}{a numeric vector of length 8 to initialize the 4 pairs of shape
and mean  for the 4 gamma distributions of the model -- 2 for burst or non-
burst durations and 2 for burst-ISI or non-burst-ISI. If specified, the fifth
and the seventh elements should match dur.shape.}
\item{hpar}{a numeric vector of length 16 specifying the 8 pairs of mean and 
sd for the 8 log-normal distributions used on the 8 gamma parameters. If 
specified, the 13th and the 15th elements should match logarithms of 
dur.shape.}
\item{tune}{a numeric vector of length 8 specifying the jump scale for the 
Metropolis proposals.}
\item{nburn}{number of iterations of the MCMC to be discarded as burn in. The 
states and parameter configurations from these iterations are not saved.}
\item{nsamp}{number of iterations of the MCMC to be saved after the burn in}
\item{nskip}{number of iterations of the MCMC to be skipped between two 
consecutive saves.}
\item{nshot}{number of Metropolis proposals to be made per update. If larger 
than 1, then corresponds to the Multuple-Try Metropolis of Liu et al (2000).}
\item{obj}{a hsmm object created by running \code{fit.hsmm}.}
\item{subset}{a subset of MCMC iterations to be used for posterior 
approximations.}
\item{glen}{grid length over which fitted ISI densities will be displayed.}
\item{newfig}{logical denoting whether plotting should be done a new x11 
device.}
\item{toplot}{logical denoting whether a plot should be generated.}
\item{hbreaks}{break points for ISI histogram.}
}
\value{
\code{fit.hsmm} returns a \var{hsmm} class object to be used by \code{print}, \code{plot} and \code{summary}.
}
\seealso{ 
\code{\link{simu.hsmm}}. 
} 
\details{Bursting and non-bursting durations are modeled as a switching gamma process with respective distributions Gamma(dur.burst.shape, dur.burst.mean) and Gamma(dur.nonburst.shape, dur.nonburst.mean). ISIs within bursts are modeled as Gamma(isi.burst.shape, isi.burst.mean) variates and similarly for non-burst ISIs. The model parameter is written as (isi.burst.shape, isi.burst.mean, isi.nonburst.shape, isi.nonburst.mean, dur.burst.shape, dur.burst.mean, dur.nonburst.shape, dur.nonburst.mean) even though (burst.shape, nonburst.shape) is treated as fixed (set equal to {dur.shape}). All these parameters are modeled as log-normal variates, entailing a total of 16 hyperparameters contained in {hpar}. Again, the hyperparameters for dur.burst.shape and dur.nonburst.shape are superfluous, but I have been too lazy to take care of these redundancies.}  
\references{Tokdar, S. T., Xi, P., Kelly, R. C. and Kass, R. E. (2009). Detection of bursts in extracellular spike trains using hidden semi-Markov point process models. \emph{Journal of Computational Neuroscience.}  DOI: 10.1007/s10827-009-0182-2

Liu, J. S., Liang, F. and Wong. W. H. (2000). The use of multiple-try 
method and local optimization in Metropolis sampling. \emph{Journal of the American Statistical Association}, 96(454): 561-573.
}

\examples{ 
## Generate a spike over 10s
## burst mean / SD in sec: (Burst) 0.2 / 0.04; (Non-burst) 0.7 / 0.07
## ISI mean / SD in sec: (Burst) 0.02 / 0.01; (Non-burst) 0.10 / 0.10
spk <- simu.hsmm(10, pars = c(25, 0.008, 100, 0.007, 4, 0.005, 1, .1))
y <- spk$spikes 

## Fit HMM
o.hmm <- fit.hsmm(y, nburn = 5e2, nsamp = 100, nskip = 5) 

## Print overview 
print(o.hmm) ## basic info about the MCMC run

## Summary of MCMC & Posterior inference
summary(o.hmm) ## nothing is printed, only trace plots of MCMC are created
sm.hmm <- summary(o.hmm, toplot = FALSE) ## will return posterior approximations

## Plot
plot(o.hmm)

## Fit HSMM
o.hsmm <- fit.hsmm(y, nburn = 5e2, nsamp = 100, nskip = 5, dur.shape = c(10,10)) 
sm.hsmm <- summary(o.hsmm, toplot = FALSE) ## will return posterior approximations
plot(o.hsmm)

## More comparative plotting
par(mfrow = c(1,2))
plot(cumsum(y), sm.hmm$prob.burst, col = 2, ty = "h", xlab = "Time stamps", ylab =  "P(Burst)", main = "HMM", bty = "n")
abline(v = cumsum(spk$transitions))
plot(cumsum(y), sm.hsmm$prob.burst, col = 2, ty = "h", xlab = "Time stamps", ylab = "P(Burst)", main = "HSMM")
abline(v = cumsum(spk$transitions), bty = "n")
} 
\keyword{programming}  