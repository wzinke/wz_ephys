library(beanplot)
library('KernSmooth')
library('MASS')

# bhvdir = '/SCRATCH/BHV_files'
# 
# bhvfl  = 'I_34_change_release_2014_06_02.dat'
# bhvfl  = 'I_34_chng_rel_2014_06_03.dat'
# bhvfl  = 'I_34_chng_rel_2014_06_04.dat'

RTbw = 25


fname  = file.choose()
bhvfl  = basename(fname)
bhvdir =  dirname(fname)

setwd(bhvdir)

dt=read.table(bhvfl, header=TRUE)

sesstime = (dt$GoEff - dt$GoEff[1])/60000

reltime   = dt$TPullEff - dt$TrialStart
relpos   = is.finite(reltime);
hitpos   = dt$Result == 0 & is.finite(dt$RT)

RThit    = density(dt$RT[hitpos], bw=RTbw, from=0, cut=0, na.rm=TRUE)
RThit$y  = RThit$y  * sum(hitpos)  * RTbw

x11(width=4, height=8, pointsize=10)

par(mfrow=c(2,1))

### reaction time distribution
max_X = tail(RThit$x, n=1)
max_Y = max(RThit$y)

truehist(dt$RT, h=RTbw, col='darkgrey', main='Reaction Time', xlab='RT [ms]', ylab='count', prob=FALSE, xlim=c(0, max_X), xaxs='i', yaxs='i')
lines(RThit,  lwd=2.5, col='red')
abline(v=median(dt$RT, na.rm=TRUE), col='red', lwd=2.5)

text(max_X, max_Y,pos=2,adj=0,labels=paste('N(hit): ',sum(hitpos), sep=''),cex=1.2)

### session time course of reaction times 
plot(sesstime[hitpos], dt$RT[hitpos], type='p', pch=20, ylab='RT [ms]', xlab='session time [min]', main='reaction times', xaxs='i', yaxs='i')
abline(h=median(dt$RT, na.rm=TRUE), col='red', lwd=2.,lty=2)

### save figure as image file
dev.copy(png,filename=paste(dt$Subject[1], '_',dt$Experiment[1],"_", dt$Date[1],'.png', sep=''),height=800, width=400, bg="white")
dev.off()



