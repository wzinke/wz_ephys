# define function to plot distribution of spike widths
MixMod = function(spkwd, histbw=25, densebw=25, main='Spike Width', x=NULL){
  require('MASS')
  require('mixtools')
  require('diptest')

  spkwd = spkwd[is.finite(spkwd)]
  N = length(spkwd)
  if(is.null(x)) {x = range(spkwd)}
  x = x[1] : x[2]

  mixmdl = normalmixEM(spkwd, k=2, maxit=5000, mu=quantile(spkwd,c(0.25,0.75)))
  print(summary(mixmdl))

  norm_1 = dnorm(x,  mixmdl$mu[1], mixmdl$sigma[1]) * mixmdl$lambda[1] * N * histbw
  norm_2 = dnorm(x,  mixmdl$mu[2], mixmdl$sigma[2]) * mixmdl$lambda[2] * N * histbw

  f <- function(val) { dnorm(val,  mixmdl$mu[1], mixmdl$sigma[1]) * mixmdl$lambda[1] -
                       dnorm(val,  mixmdl$mu[2], mixmdl$sigma[2]) * mixmdl$lambda[2] }

  distctr = sort(c(mixmdl$mu[1],mixmdl$mu[2]))

  itrsec = uniroot(f, interval=distctr, lower = distctr[1], upper=distctr[2])

  dest = density(spkwd, bw=densebw)
  dest$cnt = dest$y * N * histbw

  Dip = dest$x[dest$y == min(dest$y[dest$x > distctr[1] & dest$x < distctr[2]])]
  Dip = Dip[1]

  dipstat = dip.test(spkwd, simulate.p.value = FALSE, B = 5000)

  truehist(spkwd,h=histbw, col='grey', prob=FALSE, main=main,
           xlab=expression(paste("spike width [", mu, "s]")), ylab='count')

  lines(dest$x, dest$cnt)
  #lines(dest$x, dest$y)

  lines(x,norm_1,lwd=2,col='red')
  lines(rep(mixmdl$mu[1],2), c(0,max(norm_1)),lwd=2,col='red',lty=2)
  lines(x,norm_2,lwd=2,col='blue')
  lines(rep(mixmdl$mu[2],2), c(0,max(norm_2)),lwd=2,col='blue',lty=2)
  #lines(x,norm_1+norm_2,lwd=2,col='black',lty=3)

  abline(v=itrsec$root,lwd=2,lty=3,col='green')
  abline(v=Dip,lwd=2,lty=2,col='green')

  out = list(N=N, intersect = itrsec$root, dip=Dip, dipstat=dipstat,
             mu1=mixmdl$mu[1], sd1=mixmdl$sigma[1], lambda1=mixmdl$lambda[1],
             mu2=mixmdl$mu[2], sd2=mixmdl$sigma[2], lambda2=mixmdl$lambda[2])
}