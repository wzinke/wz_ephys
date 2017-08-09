
# compare two cumulative distributions
dist.comp = function(X1, X2, xlab=NULL, ylab='proportion', main=NULL, xlim=NULL){

  plot(range(c(X1,X2),na.rm=TRUE), c(0,1), type='n', col='red', lwd=2.5,
       xaxs='i', yaxs='i', main=main, ylab = ylab, xlab = xlab, xlim = xlim)

    sf1 = ecdf(X1)
    sf2 = ecdf(X2)
    sf3 = ecdf(X3)
    lines(sf1, col='black', do.points=FALSE, verticals=TRUE, lwd = 2)
    lines(sf2, col='red',   do.points=FALSE, verticals=TRUE, lwd = 2)
    lines(sf3, col='blue',  do.points=FALSE, verticals=TRUE, lwd = 2)


#   lines(cum.hist(X1[is.finite(X1)]), col='black', lwd=2.5)
#   lines(cum.hist(X2[is.finite(X2)]), col='red',   lwd=2.5)
#   lines(cum.hist(X3[is.finite(X3)]), col='blue',  lwd=2.5)
  lines(rep(median(X1[is.finite(X1)]),2), c(0, 0.5), col='red')
  lines(rep(median(X2[is.finite(X2)]),2), c(0, 0.5), col='blue')
  abline(h=0.5,lty=3)
  ks = ks.test(X1,X2)
  mtext(paste('p = ',ks$p.value,sep=''), side=3,adj=0)
}
