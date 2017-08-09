cum.hist = function(X){
  X = X[is.finite(X)]
  tbl = table(X, exclude=c(NA, NaN))
  ch.out   = c()
  ch.out$x = as.numeric(names(tbl))
  cnts     = as.numeric(tbl)
  ch.out$y = cumsum(cnts) / sum(cnts)

  return(ch.out)
}