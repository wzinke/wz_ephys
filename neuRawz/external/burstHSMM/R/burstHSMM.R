fit.hsmm <- function(y, scan = c(5, 25), dur.shape = c(1, 1), limit = NULL,
										 states = NULL, pars = NULL, hpar = NULL, tune = NULL,  
										 nburn = 0, nsamp = 1e3, nskip = 1, nshot = 5){
	
	
	if(length(dur.shape) == 1){
    dur.shape <- rep(dur.shape, 2)
  }
	
	n <- length(y)
  dims <- c(nspikes = n, scan.min = scan[1], scan.max = scan[2])
  mcmc_len <- c(nburn = nburn, nsamp = nsamp, nskip = nskip, nshot = nshot);

  if(is.null(hpar)){
    hpar <- c(log(10),  1,
              log(20),  2.0,
              log(10),  1,
              log(20),  2.0,
              log(dur.shape[1]),   0.5,
              log(100), 4.0,
              log(dur.shape[2]),   0.5,
              log(100), 4.0)
	}
	
  if(is.null(tune)){
    tune <- 0.15 * hpar[seq(2, 16, 2)]
	}
	
  if(is.null(pars)){
    pars <- hpar[seq(1, 16, 2)]
	}
	
  if(is.null(limit)){
    limit <- max(y) + 1
	}
  
  if(is.null(states)){
    out1 <- .C("rstate", 
							 as.integer(n), as.double(y), as.double(pars),
               states = integer(n), as.double(limit))
    
    states <- out1$states
  }


  rtime <- system.time(out <- .C("hsmm", 
																 as.double(y), as.double(pars), as.double(hpar),
                                 as.integer(states), as.integer(dims), as.integer(mcmc_len),
                                 runs = double(n), states = integer(nsamp * n),
                                 pars = double(nsamp * 8), acpt = double(4), as.double(tune),
                                 as.double(limit))
                       )[3]
  
  
  out <- list(states = matrix(out$states, nrow = n),
              runs = out$runs, 
              pars = matrix(out$pars, nrow = 8),
              runtime = rtime,
              dimensions = dims,
              mcmc = mcmc_len,
              acpt = out$acpt,
              tune = tune,
              isi = y,
              method = "HSMM"
              )
  if(any(dur.shape == 1)){
    out$method <- "HMM"
	}
    
  class(out) <- "hsmm"
  return(out)
  
}

print.hsmm <- function(obj){
  cat("Model: ", obj$method, "\n")
  print(obj$dimensions)
  print(obj$mcmc)
  cat("M-H acceptance rate = ", obj$acpt, "\n")
  cat("Run time = ", round(obj$runtime, 2),"seconds\n")
}

plot.hsmm <- function(out, subset = NULL, glen = 1e3, newfig = FALSE, hbreaks = NULL){

  if(is.null(subset)){
    subset <- 1:out$mcmc[2]
  }
	
  y <- out$isi
	if(newfig){
    x11(width = 10, height = 4)
	}
  par(mfrow = c(1, 2))

  p.states <- 1 - apply(out$states[, subset], 1, mean)
  pars <- out$pars[, subset]
  ix <- 1:8
  
  if(median(pars[2,] - pars[4,]) > 0){
    ix <- c(3, 4, 1, 2, 7, 8, 5, 6)
    p.states <- 1 - p.states
  }
  pars <- pars[ix, ]

	time.stamps <- cumsum(y)
	n.stamps <- length(y)

	burst.step <- stepfun(time.stamps, c(p.states, p.states[length(p.states)]))
	plot(burst.step, do.points = FALSE, ann = FALSE, ylim = c(-0.2, 1), axes = FALSE)
	points(time.stamps, rep(-0.1, n.stamps), pch = "|")
	axis(1, pretty(time.stamps))
	axis(2, seq(0, 1, 0.2))
  title(xlab = "Time stamps", ylab = "Burst Probability")
  g <- seq(0, 1.025 * max(y), len = glen)
  p.b <- mean(p.states)
  get.d <- function(theta, g)
    return(dgamma(g, exp(theta[1]), exp(theta[1] - theta[2])))
  
  d.b  <- p.b * apply(apply(pars[1:2, ], 2, get.d, g = g), 1, mean)
  d.nb <- (1 - p.b) * apply(apply(pars[3:4, ], 2, get.d, g = g), 1, mean)
  d.mix <- d.b + d.nb

  if(is.null(hbreaks)){
    hbreaks <- ceiling(length(y) / 5)
	}
  h <- hist(y, hbreaks, plot = FALSE)
  ymax <- max(max(h$density), max(d.mix[-1]))  
  hist(y, hbreaks, ylim = c(0, ymax), freq = FALSE,
       xlab = "ISI", main = "", col = grey(0.7))
	lines(g, d.mix, lwd = 4)
  lines(g, d.b, col = rgb(0, 1, 0))
  lines(g, d.nb, col = rgb(1, 0, 0))
	legend("top", c("Mixture", "Burst", "Non-burst"),
				 lwd = c(2, 1, 1), col = c(1, 3, 2), bty = "n")
	
  title(main = paste("Fit from", out$method, "model"), outer = TRUE, line = -2)

  invisible(cbind(grid = g, density = d.mix))
}



summary.hsmm <- function(out, subset = NULL, newfig = FALSE, toplot = TRUE){

	pars <- out$pars
	cat("dim(pars) = ", dim(pars), "\n")
  
  if(is.null(subset)){
		subset <- 1:out$mcmc[2]
  }
	
	y <- out$isi
	
  p.states <- 1 - apply(out$states[, subset], 1, mean)
  parsCI <- apply(pars[, subset], 1, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))
  ix <- 1:8
  
	if(median(pars[2,] - pars[4,]) > 0){
		ix <- c(3, 4, 1, 2, 7, 8, 5, 6)
		p.states <- 1 - p.states
	}
  parsCI <- parsCI[, ix]
	pars <- pars[ix, ]
  dimnames(parsCI)[[2]] <- c("B.ISI.shape", "B.ISI.mean", "NB.ISI.shape", "NB.ISI.mean",
														 "B.Dur.shape", "B.Dur.mean", "NB.Dur.shape", "NB.Dur.mean")
	
  if(toplot){
    if(newfig){
			x11()
		}
		par(mfrow = c(2,4))
		for(i in 1:8){
			plot(pars[i, ], ty = "l", ylab = dimnames(parsCI)[[2]][i], xlab = "Sweeps")
		}
	}
  
  invisible(list(prob.burst = p.states, par.CredInt = exp(parsCI)))
}

simu.hsmm <- function(total.time, pars = NULL){
	
  if(is.null(pars)){
    pars <- c(10, 0.01, 10, 0.05, 20, 1 / 1200, 20, 1 / 200)
	}
	
	M <- 4 * ceiling(total.time / (pars[1] * pars[2] + pars[3] * pars[4]))
	N <- 4 * ceiling(total.time / (pars[5] * pars[6] + pars[7] * pars[8]))
	
	out <- .C("rspike",
						as.double(pars), as.double(total.time), states = integer(N), spikes = double(N), 
						nspikes = integer(1), transitions = double(M), ntrans = integer(1)
						)
	
	return(list(spikes = 1e3 * out$spikes[1:out$nspikes],
							states = out$states[1:out$nspikes],
							transitions = 1e3 * out$transitions[1:out$ntrans]))
	
}

