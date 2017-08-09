\name{compModels}
\alias{compModels}
\title{Compare Duration Models on a Specific Data Set}
\description{
  Fit duration models with the maximum likelihood method to a given
  duration data set. The data can be censored. The models should be
  among the following list: inverse Gaussian, log normal, log logistic,
  gamma, Weibull, refractory exponential. The Akaike information
  criterion (AIC) is used to produce a numerical output. Diagnostic QQ
  or survival plots can also be generated.
}
\usage{
compModels(yi, ni = numeric(length(yi)) + 1,
           si = numeric(length(yi)) + 1,
           models = c("invgauss","lnorm","gamma","weibull","llogis","rexp"),
           type = c("qq","s"), log = TRUE, plot = TRUE)
}
\arguments{
  \item{yi}{vector of (possibly binned) observations or a
    \code{spikeTrain} object.}
  \item{ni}{vector of counts for each value of \code{yi}; default: \code{numeric(length(yi))+1}.}
  \item{si}{vector of counts of \emph{uncensored} observations for each
    value of \code{yi}; default: \code{numeric(length(yi))+1}.}
  \item{models}{a character vector whose elements are selected among:
    \code{"invgauss"}, \code{"lnorm"}, \code{"gamma"},
    \code{"weibull"}, \code{"llogis"}, \code{"rexp"}.}
  \item{type}{should a QQ plot (\code{"qq"}) or a survival plot
    (\code{"s"}) be generated?}
  \item{log}{should a log scale be used?}
  \item{plot}{should a plot be generated?}
}
\details{
  Fits are performed by maximizing the likelihood.
}
\value{
  A vector whose component are nammed according to the model used and
  ordered along increasing AIC values.

  if argument \code{plot} is set to \code{TRUE} (the default), a plot is
  generated as a side effect.
}
\references{
  Lindsey, J.K. (2004) \emph{The Statistical Analysis of Stochastic
    Processes in Time}. CUP.
}
\author{Christophe Pouzat  \email{christophe.pouzat@gmail.com} }
\seealso{
  \code{\link{qqDuration}},
  \code{\link{invgaussMLE}},
  \code{\link{lnormMLE}},
  \code{\link{llogisMLE}},
  \code{\link{rexpMLE}},
  \code{\link{gammaMLE}},
  \code{\link{weibullMLE}}
}
\examples{
\dontrun{
## load spontaneous data of 4 putative projection neurons
## simultaneously recorded from the cockroach (Periplaneta
## americana) antennal lobe
data(CAL1S)
## convert data into spikeTrain objects
CAL1S <- lapply(CAL1S,as.spikeTrain)
## look at the individual trains
## first the "raw" data
CAL1S[["neuron 1"]]
## next some summary information
summary(CAL1S[["neuron 1"]])
## next the renewal tests
renewalTestPlot(CAL1S[["neuron 1"]])
## It does not look too bad so let fit simple models
compModels(CAL1S[["neuron 1"]])

## Simulate a sample with 100 events from an inverse Gaussian
set.seed(1102006,"Mersenne-Twister")
mu.true <- 0.075
sigma2.true <- 3
sampleSize <- 100
sampIG <- rinvgauss(sampleSize,mu=mu.true,sigma2=sigma2.true)

## Compare models and display QQ plot
compModels(sampIG,type="qq")

## Compare models and display survival plot
compModels(sampIG,type="s")


## Generate a censored sample using an exponential distribution
sampEXP <- rexp(sampleSize,1/(2*mu.true))
sampIGtime <- pmin(sampIG,sampEXP)
sampIGstatus <- as.numeric(sampIG <= sampEXP)
## Compare models and display QQ plot
## WARNING with censored data like here the QQ plot is misleading
compModels(yi=sampIGtime,si=sampIGstatus,type="qq")
## Compare models and display survival plot
compModels(yi=sampIGtime,si=sampIGstatus,type="s")
}
}
\keyword{distribution}
\keyword{ts}
