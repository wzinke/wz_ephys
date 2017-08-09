spk_plot_histogram = function(spktm, timwin = NA, maxY=NA, max_rast=NA, EvTm=NA, binwdth=10, sm=FALSE, spikeline=FALSE, rastratio=0.25){

    require(MASS)
 
    if(is.na(timwin[1]) == 0) { spktm[spktm < timwin[1] | spktm > timwin[2]] = NA }else{ timwin = range(spktm, na.rm=TRUE) }

    ntrial = dim(spktm)[1]
    if(is.na(max_rast)) { max_rast = ntrial }

    binpos = seq(from=timwin[1]-binwdth/2,by=binwdth,to=timwin[2]+binwdth/2)
    binpos = seq(from=timwin[1],by=binwdth,to=timwin[2])

    psth = hist(spktm,breaks=binpos,plot=FALSE)

    psth$counts = (1000 * psth$counts / binwdth) / ntrial 
    maxact = ceiling(1.05*max(psth$counts))
    if(maxact > 50){
        actck = seq(0,10*floor(maxact/10),by=10*round(maxact/50))
    }else if(maxact > 25){
        actck = seq(0,5*floor(maxact/5),by=5*round(maxact/25))
    }else if(maxact > 10){
        actck = seq(0,2*floor(maxact/2),by=2*round(maxact/10))
    }else{ actck = seq(0,5*floor(maxact/5),by=round(maxact/5)) }

    plot(psth, col='black', ylim=c(0,(1+rastratio)*maxact), xlim=timwin, xaxs='i', yaxs='i', axes = F, xlab=NA, ylab=NA, main=NA)

    axis(side = 1, lwd = 2, line = -1, tck = -.02, cex=1.5, pos=0)
    axis(side = 2, lwd = 2, line = -1, las = 1, at=actck, tck = -.02, cex=1.5, pos=timwin[1])
    mtext(side = 1, 'time [ms]', line = 2)
    mtext(side = 2, 'spikes/s', line = 2.5)

    if(rastratio > 0){
        Tcnt = maxact
        tstep = (rastratio*maxact)/max_rast
        #abline(h=maxact,lwd=0.5)
        for(t in 1:ntrial) {
            Tcnt = Tcnt + tstep
            cspks = spktm[t,is.finite(spktm[t,])]
            points(cspks,rep(Tcnt, length(cspks)),pch = 20,cex=0.25)
        }
   }
}

