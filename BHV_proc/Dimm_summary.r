require(beanplot)
require('KernSmooth')
require('MASS')

# require(simecol)
# require(retimes)
# require(modeest)

RTbw = 25

# bhvdir = '/SCRATCH/BHV_files'
# 
# bhvfl = 'I_34_change_release_2014_06_02.dat'
# bhvfl = 'I_34_chng_rel_2014_06_03.dat'
# bhvfl = 'I_34_chng_rel_2014_06_04.dat'
# bhvfl = 'I_34_dimming_2014_06_05.dat'
# bhvfl = 'I_34_dimming_2014_06_06.dat'

fname  = file.choose()
bhvfl  = basename(fname)
bhvdir =  dirname(fname)

setwd(bhvdir)

dt=read.table(bhvfl, header=TRUE)

chngtime =  dt$GoEff    - dt$NoGoEff
reltime  =  dt$TPullEff - dt$NoGoEff
sesstime = (dt$NoGoEff - dt$NoGoEff[1])/60000

relpos   = is.finite(reltime);
hitpos   = dt$Result == 0 & relpos   == 1;
earlypos = dt$Result == 5 & relpos   == 1;
latepos  = dt$Result == 2 & relpos   == 1;
resppos  = hitpos    == 1 | earlypos == 1 | latepos == 1


resnm = rep(NA, length(resppos))
resnm[hitpos]   = 'accept'
resnm[earlypos] = 'early'
resnm[latepos]  = 'late'

relall   = density(reltime[resppos],         bw=RTbw, from=0, cut=0, na.rm=TRUE)
relhit   = density(reltime[hitpos],          bw=RTbw, from=0, cut=0, na.rm=TRUE)
relearly = density(reltime[earlypos],        bw=RTbw, from=0, cut=0, na.rm=TRUE)
rellate  = density(reltime[latepos],         bw=RTbw, from=0, cut=0, na.rm=TRUE)
RThit    = density(dt$RT[hitpos],            bw=RTbw, from=0, cut=0, na.rm=TRUE)
RTlate   = density(dt$RelTime[latepos],      bw=RTbw, from=0, cut=0, na.rm=TRUE)
ChngTime = density(dt$ChangeIntent[resppos], bw=RTbw, from=0, cut=0, na.rm=TRUE)

relall$y   = relall$y   * sum(relpos)   * RTbw
relhit$y   = relhit$y   * sum(hitpos)   * RTbw
relearly$y = relearly$y * sum(earlypos) * RTbw
rellate$y  = rellate$y  * sum(latepos)  * RTbw
ChngTime$y = ChngTime$y * sum(resppos)  * RTbw

RThit$y  = RThit$y  * sum(hitpos)  * RTbw
RTlate$y = RTlate$y * sum(latepos) * RTbw

# create conditions according to dimming time
#dimmitv = seq(min(dt$ChangeIntent, na.rm=TRUE), max(dt$ChangeIntent, na.rm=TRUE), length=6)
dimmitv = quantile(dt$ChangeIntent[resppos], probs = seq(0, 1, 0.2))
DimmCnd = rep(1, length(dt$ChangeIntent))

for(i in seq(length(dimmitv)-1)){
 pos = dt$ChangeIntent > dimmitv[i]
 DimmCnd[pos] = i;
}


### start plots
x11(width=12, height=8, pointsize=16)

layout(matrix(c(1, 3, 7, 2, 4, 7, 5, 5, 5, 6, 6, 6), 4, 3, byrow = TRUE))
par(mar = c(2.5, 2.5, 1.5, 0.5))
par(mgp = c(1.5, 0.5, 0))
par(oma = c(0, 0, 3, 0))

### distribution of release times
maxY = max(c(relhit$y, rellate$y , RThit$y, relearly$y,ChngTime$y))
maxX = max(c(relhit$x, RThit$x, relearly$x, rellate$x,ChngTime$x ))
plot(relhit, xlim=c(0, maxX), ylim = c(0, maxY), col = 'red', main = 'Release Time', xlab='trial time [ms]', lwd=2.5, cex=1.5, xaxs='i', yaxs='i')
#lines(relhit, col='red', lwd=2.5)
lines(relearly, col='green', lwd=2.5)
lines(rellate, col='blue', lwd=2.5)
lines(ChngTime, col='black', lwd=1)

abline(v=dimmitv, lty=2, col='black')

legend("topright", c('hit','early','late','change'), pch = 15, col = c('red','green','blue','black'),ncol=1,bty='n',cex=1.2)


### RT distribution
max_X = max(c(tail(RThit$x, n=1), tail(RTlate$x, n=1)))
max_Y = max(c(max(RThit$y), max(RTlate$y)))

truehist(dt$RT, h=RTbw, col='darkgrey', main='Reaction Time', xlab='RT [ms]', ylab='count', prob=FALSE, xlim=c(0, max_X), xaxs='i', yaxs='i')
lines(RThit,  lwd=2.5, col='red')
lines(RTlate, lwd=2.5, col='blue')
abline(v=median(dt$RT, na.rm=TRUE), col='red', lwd=2.5)

text(max_X, max_Y/2,pos=2,adj=0,labels=paste('N(hit): ',sum(hitpos),'\n','N(early): ',sum(earlypos),'\n','N(late): ',sum(latepos),'\n',sep=''),cex=1.2)

## RT for dimming intervals
# beanplot(dt$RT[hitpos]~DimmCnd[hitpos], beanlines='median', overallline='median', b w=RTbw, col='darkgrey', main='Reaction Times', xlab='dimming interval', ylab='RT [ms]', yaxs='i')
boxplot(dt$RT[hitpos]~DimmCnd[hitpos], notch=TRUE, outline=FALSE, col='darkgrey', main='Reaction Times', xlab='dimming interval', ylab='RT [ms]', yaxs='i', ylim=range(dt$RT[hitpos], na.rm=TRUE))
stripchart(dt$RT[hitpos]~DimmCnd[hitpos], method = "jitter", jitter = 0.1, vertical =TRUE, add=TRUE, pch=20, cex=0.8)
abline(h=median(dt$RT[hitpos]), lty=2, lwd=2, col='red')

### performance for dimming intervals
barplot(100 * prop.table(table(resnm[resppos], DimmCnd[resppos]), margin=2), beside=T, space=c(0, 0.4), col=c('red', 'green', 'blue'), ylab='rate [%]', xlab='dimming interval', main='Performance', yaxs='i')

### session time course of release times 
rngX = range(sesstime[resppos], na.rm=TRUE)
rngY = range(reltime[resppos], na.rm=TRUE)

plot(sesstime[hitpos], reltime[hitpos], type='p', pch=20, ylab='release time [ms]', xlab='session time [min]', main='Release times', col='red', xaxs='i', yaxs='i', xlim=rngX, ylim=rngY)

points(sesstime[earlypos], reltime[earlypos], pch=20, col='green')
points(sesstime[latepos],  reltime[latepos],  pch=20, col='blue')

### session time course of reaction times 
plot(sesstime[hitpos], dt$RT[hitpos], type='p', pch=20, ylab='RT [ms]', xlab='session time [min]', main='reaction times', xaxs='i', yaxs='i', xlim=rngX)
abline(h=median(dt$RT, na.rm=TRUE), col='red', lwd=2.,lty=2)

### release times dependent on dimming time
rngX = range(dt$ChangeIntent[resppos], na.rm=TRUE)
# rngY = range(reltime[resppos], na.rm=TRUE)

plot(dt$ChangeIntent[hitpos], reltime[hitpos], type='p', pch=20, ylab='release time [ms]', xlab='change time [ms]', main='Release times', col='red', xaxs='i', yaxs='i', xlim=rngX, ylim=rngY)
points(dt$ChangeIntent[earlypos], reltime[earlypos], pch=20, col='green')
points(dt$ChangeIntent[latepos], reltime[latepos], pch=20, col='blue')

abline(a=0, b=1, , lty=3, lwd=2)

### figure title
mtext(paste(dt$Subject[1], " - Performance Summary from ", dt$Date[1], sep=''), outer = TRUE, side = 3, cex = 1.2, line = 1)

### save figure as image file
dev.copy(png,filename=paste(dt$Subject[1], '_',dt$Experiment[1],"_", dt$Date[1],'.png', sep=''),height=800, width=1200, bg="white")
dev.off()
