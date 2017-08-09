DenseHist = function(X,main=NULL, xlab=NULL,ylab=NULL, h=NULL, xlim=NULL){
  require(MASS)

  X= X[is.finite(X)]

  if(is.null(xlim)) { xlim = range(X) }

  truehist(X,h=h, col='grey', prob=FALSE, main=main,
           xlab=xlab, ylab=ylab, xlim=xlim)

   dest = density(X, bw=h)
   dest$cnt = dest$y * length(X) * h
   lines(dest$x, dest$cnt,lwd=2,col='blue')
   mtext(paste('N = ',length(X),sep=''), side=3,adj=0)
}
