\name{crossGeneral}
\alias{crossGeneral}
\alias{crossTight}
\alias{mkTightBMtargetFct}
\alias{print.FirstPassageTime}
\alias{summary.FirstPassageTime}
\alias{plot.FirstPassageTime}
\alias{lines.FirstPassageTime}
\title{Computations of Boundary Crossing Probabilities for the Wiener Process}
\description{
  Computes the distribution of the first passage time through an arbitrary
  (\code{crossGeneral}) or a "tight" (\code{crossTight}) boundary for a Wiener
  process. The method of Loader and Deely (1987) is used. A tight
  boundary is a boundary generating the tighest confidence band for the
  process (Kendall et al, 2007). Utility function and methods:
  \code{mkTightBMtargetFct}, \code{print}, \code{summary}, \code{plot},
  \code{lines}, are also provided to use and explore the results. 
}
\usage{
crossGeneral(tMax = 1, h = 0.001, cFct, cprimeFct, bFct, withBounds = FALSE, Lplus)
crossTight(tMax = 1, h = 0.001, a = 0.3, b = 2.35, withBounds = TRUE, logScale = FALSE)
mkTightBMtargetFct(ci = 0.95, tMax = 1, h = 0.001, logScale = FALSE)
\method{print}{FirstPassageTime}(x, ...)
\method{summary}{FirstPassageTime}(object, digits, ...)
\method{plot}{FirstPassageTime}(x, y, which = c("Distribution", "density"), xlab, ylab, ...)
\method{lines}{FirstPassageTime}(x, which = c("Distribution", "density"), ...)
}
\arguments{
  \item{tMax}{A positive \code{numeric}. The "time" during which the
    Wiener process is followed.}
  \item{h}{A positive \code{numeric}. The integration time step used for
    the numerical solution of the Volterra integral equation (see
    \code{details}).
  }
  \item{cFct}{A \code{function} defining the boundary to be crossed. The
    first argument of the function should be a "time" argument. If the
    first argument is a vector, the function should return a vector of the
    same length.
  }
  \item{cprimeFct}{A \code{function} defining time derivative of the
    boundary to be crossed. Needs to be specified only if a check of
    the sign of the kernel derivative (see \code{details}) is requested. The
    first argument of the function should be a "time" argument. If the
    first argument is a vector, the function should return a vector of the
    same length.
  }
  \item{bFct}{A \code{function}. The "b" function of Loader and Deely
    (1987). Does not need to be specified (\emph{i.e.}, can be
    missing) but can be used to improve convergence. The
    first argument of the function should be a "time" argument. If the
    first argument is a vector, the function should return a vector of the
    same length.
  }
  \item{withBounds}{A \code{logical}. Should bounds on the distribution
    be calculated? If yes, set it to \code{TRUE}, leave it to its
    default value, \code{FALSE}, otherwise.}
  \item{Lplus}{A \code{logical}. If bounds are requested
    (\code{withBounds=TRUE}) and if the sign of the time derivative of
    the kernel is known to be positive or null, set to \code{TRUE}, if
    it is known to be negative, set it to \code{FALSE}. If the sign is
    unknown, leave \code{Lplus} unspecified and provide a
    \code{cprimeFct} function.
  }
  \item{logScale}{A \code{logical}. Should intermediate calculations in
    \code{crossTight} be carried out on the log scale for numerical
    precision? If yes, set it to \code{TRUE}, leave it to its default,
    \code{FALSE}, otherwise.}
  \item{a,b}{\code{numerics}, the two parameters of the "tight"
    boundary: \code{c(t) = a + b*sqrt(t)}. See details.}
  \item{ci}{A \code{numeric} larger than 0 and smaller than 1. The
    nominal coverage probability desired for a "tight" confidence band
    (see \code{details}).
  }
  \item{x,object}{A \code{FirstPassageTime} object returned by
    \code{crossGeneral} or \code{crossTight}.
  }
  \item{y}{Not used but required for a \code{plot} method.}
  \item{which}{A \code{character} string, "\code{Distribution}" or
    "\code{density}", specifying if a probability distribution or a
    probability density should be graphed.
  }
  \item{xlab,ylab}{See \code{\link{plot}}.}
  \item{digits}{A positive \code{integer}. The number of digits to print
    in \code{summary}. If bounds were computed, the value of \code{digits}
    is computed internally based on the bounds width.
  }
  \item{\dots}{Used in \code{plot} and \code{lines} to pass further
    arguments (see \code{\link{plot}} and \code{\link{lines}}), not used
    in \code{print} and \code{summary}.
  }
  
}
\details{
  The Loader and Deely (1987) method to compute the probability
  \eqn{G(t)}{G(t)} that the
  first passage of a Wiener process / Brownian motion occurs between 0
  and \eqn{t}{t} (argument \code{tMax} of \code{crossGeneral} and
  \code{crossTight}) through a boundary defined by \eqn{c(t)}{c(t)} is based
  on the numerical solution of a Volterra 
  integral equation of the first kind satisfied by \eqn{G()}{G()} and
  defined by their Eq. 2.2:
  \deqn{F(t) = \int_0^t K(t,u) dG(u)}{F(t) = \int_0^t K(t,u) dG(u)}
  where, \eqn{F(t)}{F(t)} is defined by:
  \deqn{F(t)=\Phi(-\frac{c(t)}{\sqrt{t}})+\exp \big( -2 b(t) \,
  (c(t)-tb(t))\big) \, \Phi(\frac{-c(t)+2\, t \, b(t)}{\sqrt{t}})}{F(t)
  = pnorm(-c(t)/sqrt(t)) +
  exp(-2*b(t)*(c(t)-t*b(t)))*pnorm((-c(t)+2*t*b(t))/sqrt(t))}
  \eqn{K(t,u)}{K(t,u)} is defined by:
  \deqn{K(t,u)=\Phi(\frac{c(u)-c(t)}{\sqrt{t-u}})+\exp \big(-2 b(t) \, (c(t)
    -c(u) -(t-u) b(t))\big) \, \Phi(\frac{c(u)-c(t)+2\, (t-u) \,
    b(t)}{\sqrt{t-u}})}{K(t,u)=pnorm((c(u)-c(t))/sqrt(t-u)) +
    exp(-2*b(t)*(c(t)-c(u)-(t-u)*b(t)))*pnorm((c(u)-c(t)+2*(t-u)*b(t))/(sqrt(t-u)))
  }
  and \eqn{b(t)}{b(t)} is an additional function (that can be uniformly
  0) that is chosen to improve convergence speed and to get error
  bounds. Argument \code{h} is the step size used in the numerical
  solution of the above Volterra integral equation. The mid-point method
  (Eq. 3.1 and 3.2 of Loader and Deely (1987)) is implemented. If
  \code{tMax} is not a multiple of \code{h} it is modified as follows:
  \code{tMax <- round(tMax/h)*h}.

  \code{crossGeneral} generates functions \eqn{F()}{F()} and \eqn{K(,)}{K(,)}
  internally given \eqn{c()}{c()} (argument \code{cFct}) and \eqn{b()}{b()}
  (argument \code{bFct}). If \code{bFct} is not given (\emph{i.e.},
  \code{missing(bFct)} returns \code{TRUE}) it is taken as uniformly
  0. If a numeric is given for \code{cFct} then \eqn{c()}{c()} is
  defined as a uniform function returning the first element of the
  argument (\code{cFct}).

  Function \code{crossTight} assumes the following functional form for
  \eqn{c()}{c()}: \eqn{c(t) = a + b \, \sqrt(t)}{a + b *
  sqrt(t)}. \eqn{b()}{b()} is set to \eqn{c'()}{c'()} (the derivative of
  \eqn{c()}{c()}). Arguments \code{a} and \code{b} of \code{crossTight}
  correspond to the 2 parameters of \eqn{c()}{c()}.

  If argument \code{withBounds} is set to \code{TRUE} then bounds on
  \eqn{G()}{G()} are computed. Function \code{crossTight} uses Eq. 3.6
  and 3.7 of Loader and Deely (1987) to compute these bounds,
  \eqn{G_{L}(t)}{Gl(t)} and \eqn{G_U(t)}{Gu(t)}. Function
  \code{crossGeneral} uses Eq. 3.6 and 3.7 (if argument \code{Lplus} is
  set to \code{TRUE}) or Eq. 3.10 and 3.11 (if argument \code{Lplus} is
  set to \code{FALSE}). Here \code{Lplus} stands for the sign of the
  partial derivative of the kernel \eqn{K(,)}{K(,)} with respect to its
  second argument. If the sign is not known the user can provide the
  derivative \eqn{c'()}{c'()} of \eqn{c()}{c()} as argument
  \code{cprimeFct}. A (slow) numerical check is then performed to decide
  wether \code{Lplus} should be \code{TRUE} or \code{FALSE} or if it
  changes sign (in which case bounds cannot be obtained and an error is returned).

  In function \code{crossTight} argument \code{logScale} controls the
  way some intermediate computations of the mid-point method are
  implemented. If set to \code{FALSE} (the default) a literal
  implementation of Eq. 3.2 of Loader and Deely (1987) is used. If set
  to \code{TRUE} then additions subtractions are computed on the log
  scale using functions \code{logspace_add} and \code{logspace_sub} of
  the \code{R API}. The computation is then slightly slower and it turns
  out that the gain
  in numerical precision is not really significant, so you can safely
  leave this argument to its default value.
}
\value{
  \code{crossGeneral} and \code{crossTight} return a
  \code{FirstPassageTime} object which is a \code{list} with the
  following components:
  \item{time }{A \code{numeric} vector of "times" at which the first
    passage time probability has been evaluated.}
  \item{G }{A \code{numeric} vector of first passage probability.}
  \item{Gu }{A \code{numeric} vector with the upper bound of first
    passage probability. Only if \code{withBounds} was set to
    \code{TRUE}.
  }
  \item{Gl }{A \code{numeric} vector with the lower bound of first
    passage probability. Only if \code{withBounds} was set to
    \code{TRUE}.
  }
  \item{mids }{A \code{numeric} vector of "times" at which the first
    passage time probability \emph{density} has been evaluated. Mid
    points of component \code{time}.}
  \item{g }{A \code{numeric} vector of first passage probability
    \emph{density}.}
  \item{h }{A \code{numeric}. The value of argument \code{h} of
    \code{crossGeneral} or \code{crossTight}.}
  \item{call }{The matched call.}

  \code{mkTightBMtargetFct} returns a \code{function} which can be used
  in optim. This function returns the square of the difference between
  \code{(1-ci)/2} (remember the "symmetry" of the Wiener processes
  paths, that is, for every path there is a symmetric one with respect
  to the abscissa having \emph{with the same probability}) and the
  probability to have the first passage time of 
  the Wiener process smaller or equal to 1 when the boundary is the
  "tight" boundary defined by: \eqn{a + b \, \sqrt{t}}{a +
  b*sqrt(t)}. The function takes a single vector argument containing the
  \emph{log} of the parameters \code{a} (vector's first element) and
  \code{b} (vector's second element).
  
  
  Methods \code{print.FirstPassageTime} and
  \code{summary.FirstPassageTime} output the probability to observe the
  first exit between 0 and \code{tMax}. If bounds were computed, the
  precision on the probability is also returned (as an attribute for
  \code{print.FirstPassageTime}). \code{summary.FirstPassageTime} also
  gives the integration time step, \code{h}, used.
}
\references{
  C. R. Loader and J. J. Deely (1987) Computations of Boundary
    Crossing Probabilities for the Wiener
    Process. \emph{J. Statist. Comput. Simul.} \bold{27}: 95--105.

  W. S. Kendall, J.- M. Marin and C. P. Robert (2007) Brownian
  Confidence Bands on Monte Carlo Output. \emph{Statistics and
  Computing} \bold{17}: 1--10. Preprint available at: \url{http://www.ceremade.dauphine.fr/\%7Exian/kmr04.rev.pdf}
}
\author{ Christophe Pouzat  \email{christophe.pouzat@gmail.com} }
\note{
  Using \code{logScale = TRUE} in \code{crossTight} seems to be an
  overkill, that is, it doubles computation's time without bringing
  significant numerical improvement.

  \code{crossGeneral} is for now pure \code{R} code. The first passage
 probability is obtained by solving the lower triangular system (Eq. 3.1
 of Loader and Deely (1987)) with \code{\link{forwardsolve}} and is
 therefore rather fast (but can be memory "hungry"). The bounds are
 computed by 2 nested loops and can therefore be long to compute.

 \code{crossTight} is calling a \code{C} code and is fast.

 Loader and Deely paper also describes a method where \eqn{G(t)}{G(t)}
 is solution of a Volterra integral equation of the second kind (their
 Eq. 2.7). This approach is presently not implemented in the above
 functions.
 
}
\section{Warning }{
  \code{crossGeneral} with \code{withBounds = TRUE} and a negative
  kernel derivative is presently poorly tested, so be careful and let me
  know if mistakes show up.
} 
\seealso{
  \code{\link{print}},
  \code{\link{summary}},
  \code{\link{plot}},
  \code{\link{lines}},
  \code{\link{pinvgauss}}
}
\examples{
\dontrun{
## Reproduce Table 1 (p 101) of Loader and Deely (1987)
## define a vector of n values
nLD <- c(8,16,32,64,128)

## Part 1: c(t) = sqrt(1+t) and tMax=1
## define cFct
cFT1p1 <- function(t) sqrt(1+t)
## define the different bFct
bFT1p1.ii <- function(t) 0.5/sqrt(1+t)
bFT1p1.iii <- function(t) (cFT1p1(t)-cFT1p1(0))/t 
bFT1p1.iv <- function(t) 0.5*(bFT1p1.ii(t)+bFT1p1.iii(t)) 
bFT1p1.v <- function(t) (2*t-4/5*((1+t)^2.5-1))/t^3+3*cFT1p1(t)/2/t
## Do the calculations
round(t(sapply(nLD,
               function(n) {
                 c(n=n,
                   i=crossGeneral(tMax=1,h=1/n,cFct=cFT1p1)$G[n],
                   ii=crossGeneral(tMax=1,h=1/n,cFct=cFT1p1,bFct=bFT1p1.ii)$G[n],
                   iii=crossGeneral(tMax=1,h=1/n,cFct=cFT1p1,bFct=bFT1p1.iii)$G[n],
                   iv=crossGeneral(tMax=1,h=1/n,cFct=cFT1p1,bFct=bFT1p1.iv)$G[n],
                   v=crossGeneral(tMax=1,h=1/n,cFct=cFT1p1,bFct=bFT1p1.v)$G[n])})),
      digits=6)

## Part 2: c(t) = exp(-t) and tMax=1
## define cFct
cFT1p2 <- function(t) exp(-t)
## define the different bFct
cFT1p2 <- function(t) exp(-t)
bFT1p2.ii <- function(t) -exp(-t)
bFT1p2.iii <- function(t) (cFT1p2(t)-cFT1p2(0))/t 
bFT1p2.iv <- function(t) 0.5*(bFT1p2.ii(t)+bFT1p2.iii(t)) 
bFT1p2.v <- function(t) 3*(1-t-exp(-t))/t^3+3*cFT1p2(t)/2/t
## Do the calculations
round(t(sapply(nLD,
               function(n) {
                 c(n=n,
                   i=crossGeneral(tMax=1,h=1/n,cFct=cFT1p2)$G[n],
                   ii=crossGeneral(tMax=1,h=1/n,cFct=cFT1p2,bFct=bFT1p2.ii)$G[n],
                   iii=crossGeneral(tMax=1,h=1/n,cFct=cFT1p2,bFct=bFT1p2.iii)$G[n],
                   iv=crossGeneral(tMax=1,h=1/n,cFct=cFT1p2,bFct=bFT1p2.iv)$G[n],
                   v=crossGeneral(tMax=1,h=1/n,cFct=cFT1p2,bFct=bFT1p2.v)$G[n])})),
      digits=6)

## Part 3: c(t) = t^2 + 3*t + 1 and tMax=1
## define cFct
cFT1p3 <- function(t) t^2+3*t+1
## define the different bFct
bFT1p3.ii <- function(t) 2*t+3
bFT1p3.iii <- function(t) (cFT1p3(t)-cFT1p3(0))/t 
bFT1p3.v <- function(t) 5*t/4+3
bFT1p3.vi <- function(t) rep(3,length(t))
round(t(sapply(nLD,
               function(n) {
                 c(n=n,
                   i=crossGeneral(tMax=1,h=1/n,cFct=cFT1p3)$G[n],
                   ii=crossGeneral(tMax=1,h=1/n,cFct=cFT1p3,bFct=bFT1p3.ii)$G[n],
                   iii=crossGeneral(tMax=1,h=1/n,cFct=cFT1p3,bFct=bFT1p3.iii)$G[n],
                   v=crossGeneral(tMax=1,h=1/n,cFct=cFT1p3,bFct=bFT1p3.v)$G[n],
                   vi=crossGeneral(tMax=1,h=1/n,cFct=cFT1p3,bFct=bFT1p3.vi)$G[n])})),
      digits=6)

## Part 3: c(t) = t^2 + 3*t + 1 and tMax=1
## define cFct
cFT1p4 <- function(t) 1/t
## Here only column (i) and (vii) are reproduced.
## If one attempts to reproduce (ii) directly with crossGeneral
## a NaN appears (when a -Inf is the correct value) in functions
## F() and K(,) (see details) for instance when t=0 in F.
## Then as crossGeneral is presently written R attempts to
## compute t*b(t), where b(t) is c'(t), that is, t*(-1/t^2) which is
## NaN (for R) when t=0.
bFT1p4.vii <- function(t) rep(-1,length(t))
round(t(sapply(nLD,
               function(n) {
                 c(n=n,
                   i=crossGeneral(tMax=1,h=1/n,cFct=cFT1p4)$G[n],
                   vii=crossGeneral(tMax=1,h=1/n,cFct=cFT1p4,bFct=bFT1p4.vii)$G[n])})),
      digits=6)
## The last 3 rows of column (vii) are not the same as in the paper

## Reproduce Table 4 (p 104) of Loader and Deely (1987)
## As before the probability of first passage between
## 0 and 1 is computed. This time only three boundary
## functions are considered. The error bounds are
## obtained

## Part 1: c(t) = sqrt(1+t)
## Left columns pair: b(t) = c'(t)
round(t(sapply(nLD,
               function(n) {
                 res <- crossGeneral(tMax=1,h=1/n,cFct=cFT1p1,bFct=bFT1p1.ii,withBounds=TRUE,Lplus=TRUE)
                 c(Gl=res$Gl[n],Gu=res$Gu[n])
               }
               )
         ),
       digits=5)

## Right columns pair: b(t) = 0
round(t(sapply(nLD,
               function(n) {
                 res <- crossGeneral(tMax=1,h=1/n,cFct=cFT1p1,withBounds=TRUE,Lplus=TRUE)
                 c(n=n,Gl=res$Gl[n],Gu=res$Gu[n])
               }
               )
         ),
       digits=5)

## Part 2: c(t) = t^2 + 3*t + 1
## Left columns pair: b(t) = 3 - 2*t
round(t(sapply(nLD,
               function(n) {
                 res <- crossGeneral(tMax=1,h=1/n,cFct=cFT1p3,bFct=function(t) 3-2*t,withBounds=TRUE,Lplus=TRUE)
                 c(n=n,Gl=res$Gl[n],Gu=res$Gu[n])
               }
               )
         ),
       digits=5)

## Right columns pair: b(t) = 3 - t
round(t(sapply(nLD,
               function(n) {
                 res <- crossGeneral(tMax=1,h=1/n,cFct=cFT1p3,bFct=function(t) 3-2*t,withBounds=TRUE,Lplus=TRUE)
                 c(n=n,Gl=res$Gl[n],Gu=res$Gu[n])
               }
               )
         ),
       digits=5)

## Part 3: c(t) = 1 + sin(t)
## Left columns pair: b(t) = c'(t)
round(t(sapply(nLD,
               function(n) {
                 res <- crossGeneral(tMax=1,h=1/n,cFct=function(t) 1+sin(t),bFct=function(t) cos(t),withBounds=TRUE,Lplus=TRUE)
                 c(n=n,Gl=res$Gl[n],Gu=res$Gu[n])
               }
               )
        ),
      digits=5)

## Left columns pair: b(t) = 0.5
round(t(sapply(nLD,
               function(n) {
                 res <- crossGeneral(tMax=1,h=1/n,cFct=function(t) 1+sin(t),bFct=function(t) rep(0.5,length(t)),withBounds=TRUE,Lplus=TRUE)
                 c(n=n,Gl=res$Gl[n],Gu=res$Gu[n])
               }
               )
        ),
      digits=5)


## Check crossGeneral against an analytical inverse Gaussian
## distribution
## Define inverse Gaussian parameters
mu.true <- 0.075
sigma2.true <- 3
## Define a function transforming the "drift" (mu.true) and
## "noise variance" (sigma2.true) of the default inverse
## Gaussian parametrization of STAR into a
## linear boundary of an equivalent Wiener process first
## passage time problem
star2ld <- function(mu,sigma2) c(sqrt(1/sigma2),-sqrt(1/sigma2)/mu)
## Get the "equivalent" boundary parameters (y intercept and slope)
parB1 <- star2ld(mu.true,sigma2.true)
## Plot the "target" inverse Gaussian density
xx <- seq(0.001,0.3,0.001)
plot(xx,dinvgauss(xx,mu=mu.true,sigma2=sigma2.true),type="l")
## Get the numerical estimate of the density using Loader and
## Deely Volterra integral equation method
igB1 <- crossGeneral(tMax=0.3,h=0.001,cFct=function(t) parB1[1]+parB1[2]*t,withBounds=FALSE)
## superpose the numerical estimate to the exact solution
## use lines method to do that
lines(igB1,"density",col=2)

## Use of crossTight and associated function
## Get the paramters, a and b, of the "approximate"
## tightest boundary: c(t) = a + b*sqrt(t), giving a 
## 0.05 probability of exit between 0 and 1
## (in fact we are discussing here a pair of symmetrical
## bounaries, c(t) and -c(t)). See Kendall et al (2007)
## for details
## Start by defining the target function
target95 <- mkTightBMtargetFct(ci=0.95)
## get the optimal log(a) and log(b) using
## the values of table 1 of Kendall et al as initial
## guesses
p95 <- optim(log(c(0.3,2.35)),target95,method="BFGS")
## check the convergence of BFGS
p95$convergence
## check if the parameters changed a lot
exp(p95$par)
## Get the bounds on G(1) for these optimal parameters
d95 <- crossTight(a=exp(p95$par[1]),b=exp(p95$par[2]),withBound=TRUE,logScale=FALSE)
## print out the summary
summary(d95)
## Do the same for the 0.01 probability of first passage
target99 <- mkTightBMtargetFct(ci=0.99)
p99 <- optim(p95$par,target99,method="BFGS")
p99$convergence
exp(p99$par)
d99 <- crossTight(a=exp(p99$par[1]),b=exp(p99$par[2]),withBound=TRUE,logScale=FALSE)
summary(d99) 
}

}
\keyword{distribution}
\keyword{htest}

