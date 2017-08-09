library(beanplot)
library('KernSmooth')
library('MASS')
#library(simecol)
library(retimes)
library(modeest)

bhvdir = '/SCRATCH/BHV_files'

bhvfl  = 'Ingrid_get_joystick_2014_04_26.dat'
bhvfl  = 'Ingrid_get_joystick_2014_04_28.dat'
bhvfl  = 'I_34_change_release_2014_06_02.dat'
bhvfl  = 'I_34_chng_rel_2014_06_03.dat'
bhvfl  = 'I_34_chng_rel_2014_06_04.dat'
bhvfl  = 'I_34_dimming_2014_06_05.dat'
bhvfl  = 'I_34_dimming_2014_06_06.dat'
# fname <- file.choose()

setwd(bhvdir)

dt=read.table(bhvfl, header=TRUE)
RTbw = 25

reltime = dt$TPullEff - dt$TrialStart
relpos   = is.finite(reltime);
hitpos   = dt$Result == 0 & relpos == 1;
earlypos = dt$Result == 5 & relpos == 1;
latepos  = dt$Result == 2 & relpos == 1;
resppos  = hitpos    == 1 | earlypos == 1 | latepos == 1


relall   = density(dt$RelTime[resppos],  bw=RTbw, from=0, cut=0, na.rm=TRUE)
relhit   = density(dt$RelTime[hitpos],  bw=RTbw, from=0, cut=0, na.rm=TRUE)
relearly = density(dt$RelTime[earlypos],bw=RTbw, from=0, cut=0, na.rm=TRUE)
RThit    = density(dt$RT[hitpos],           bw=RTbw, from=0, cut=0, na.rm=TRUE)
ChngTime = density(dt$ChangeIntent[resppos],  bw=RTbw, from=0, cut=0, na.rm=TRUE)


 relall$y   = relall$y   * sum(relpos)
 relhit$y   = relhit$y   * sum(hitpos)
 relearly$y = relearly$y * sum(earlypos)
 ChngTime$y = ChngTime$y * sum(resppos)

# RTall$y   = RTall$y   * sum(relpos)
# RThit$y   = RThit$y   * sum(hitpos)

maxY = max(c(relall$y,  relhit$y, RThit$y, relearly$y))
maxX = max(c(c(relall$x, relhit$x, RThit$x, relearly$x)))


# create conditions according to dimming time
dimmitv = seq(min(dt$DimmIntend), max(dt$DimmIntend),length=5)

dt$DimmCnd = rep(1,length(dt$DimmIntend))

for(i in seq(length(dimmitv)-1)){
    pos = dt$DimmIntend > dimmitv[i]
    dt$DimmCnd[pos] = i;
}



x11()

#par(mfrow=c(2,3))

layout(matrix(c(1,2,3,4,4,4,5,5,5), 3, 3, byrow = TRUE))




#x11()

plot(density(dt$RT,na.rm=TRUE,bw=RTbw),lwd=2.5,col='blue',main='release times',xlab='trial time [ms]')

truehist(dt$RT,h=25,col='darkgrey',main='response times',xlab='RT [ms]')
lines(density(dt$RT,na.rm=TRUE,bw=RTbw),lwd=2.5,col='blue')
abline(v=median(dt$RT,na.rm=TRUE),col='red',lwd=2.5)



plot((dt$TrialStart-dt$TrialStart[1])/60000,dt$RT,type='p',pch=20,ylab='RT [ms]',xlab='session time [min]')




dt$Tstart  = dt$Tstart  * 1000;
dt$TFixOn  = dt$TFixOn  * 1000 - dt$Tstart 
dt$WaitEff = abs(dt$TFixOn)
dt$DimmEff = dt$DimmEff * 1000 - dt$Tstart 
dt$TRel    = dt$TRel    * 1000 - dt$Tstart 

stpos    = is.finite(dt$WaitEff);
relpos   = is.finite(dt$ReleaseTime);
hitpos   = dt$Result == 0 & relpos == 1;
earlypos = dt$Result == 5 & relpos == 1;
resppos  = hitpos    == 1 | earlypos == 1

# 
# x11()
# timefit(dt$WaitEff[stpos],plot=TRUE,iter=1000)
# x11()
# 
# timefit(dt$RT[hitpos],plot=TRUE,iter=1000)


x11()

par(mfrow=c(2,3))
RTbw = 50;

mlv(dt$WaitEff[stpos], method = "mfv")

    truehist(dt$WaitEff[stpos], col='black',border='grey',h=RTbw,prob=TRUE,xlab='time [ms]',ylab='density',main='Wait Start')
#     plot(Tstart, col = 'black',main = 'Trial Start Times',xlab='trial time [ms]',lwd=2.5, cex=1.5)

    lines(Tstart, col = 'red',lwd=2.5, cex=1.5)
    rug(dt$WaitStart[stpos], side = 1, col = "black")


    plot(relall,xlim=c(0,maxX), ylim = c(0,maxY), col = 'black',main = 'Release times',xlab='trial time [ms]',lwd=2.5, cex=1.5)
    lines(relhit,col='red',lwd=2.5)
    lines(relearly,col='green',lwd=2.5)
    lines(ChngTime,col='blue',lwd=1.5)

    rug(dt$DimmTime, side = 1, col = "blue")
    # rug(jitter(relearly, amount = 0.01), side = 3, col = "green")

    truehist(dt$RT[hitpos], col='black',border='grey',h=RTbw,prob=TRUE,xlab='time [ms]',ylab='density',main='Response time')
#     plot(RThit,xlim=c(0,maxX), ylim = c(0,maxY), col = 'red',main = 'Reaction times',xlab='response time [ms]',lwd=2.5, cex=1.5)
     lines(RThit, col = 'red',lwd=2.5, cex=1.5)
   rug(dt$RT[hitpos], side = 1, col = "black")

RTbw=0.05
    beanplot(ReleaseTime ~ DimmCnd, data = dt, subset = hitpos, col = "salmon", border = "grey", cutmin = 0, bw = RTbw, overallline = "median", beanlines= "median",side='second', xlim=c(0.5,max(dt$DimmCnd)+0.5), ylim=range(dt$ReleaseTime[relpos]),  at=sort(unique(dt$DimmCnd[hitpos])),main='Release times',ylab='time [ms]') 
    beanplot(ReleaseTime ~ DimmCnd, data = dt, subset = earlypos, col = "lightgreen", border = "grey", cutmin = 0, bw = RTbw, overallline = "median", beanlines= "median",side='first', add=TRUE, at=sort(unique(dt$DimmCnd[earlypos])))

RTbw=50

    beanplot(RT ~ DimmCnd, data = dt, subset = hitpos, col = "salmon", border = "grey", cutmin = 0, bw = RTbw, overallline = "median", beanlines= "median",main='Response times',ylab='time [ms]') 

    barplot(prop.table(table(dt$Result,dt$DimmCnd),2),beside=TRUE,col=c('salmon','lightblue','lightgreen'),ylab='rate [%]',main='performance')

#     lines(relearly,col='green',lwd=2.5)
#     rug(dt$ReleaseTime[earlypos]], side = 1, col = "red")


# rug(jitter(WIboot, amount = 0.01), side = 1, col = "red")
# rug(jitter(ACboot, amount = 0.01), side = 3, col = "blue")


#     abline(v=mean(ACboot),col='blue',lty=3,lwd=2.5)
#     abline(v=mean(WIboot),col='red',lty=3,lwd=2.5)
#     abline(v=mean(DFboot),col='green',lty=3,lwd=2.5)


x11()
par(mfrow=c(4,1))




plot(dt$TrialNo[hitpos]  ,dt$ReleaseTime[hitpos], type='p',pch=20, col='red',ylim=c(range(dt$ReleaseTime,na.rm=TRUE)), xlim=c(range(dt$TrialNo,na.rm=TRUE)),xaxs='i')
lines(locpoly(dt$TrialNo[hitpos]  ,dt$ReleaseTime[hitpos],bandwidth=RTbw), col='darkred')
mod=rlm(dt$ReleaseTime[hitpos]~dt$TrialNo[hitpos] )
abline(mod,lty=3,lwd=2,col='darkred')


points(dt$TrialNo[earlypos]  ,dt$ReleaseTime[earlypos], pch=20, col='green')
lines(locpoly(dt$TrialNo[earlypos]  ,dt$ReleaseTime[earlypos],bandwidth=RTbw), col='darkgreen')
mod=rlm(dt$ReleaseTime[earlypos]~dt$TrialNo[earlypos] )
abline(mod,lty=3,lwd=2,col='darkgreen')




plot(dt$TrialNo[hitpos]  ,dt$RT[hitpos], type='p',pch=20, col='red',ylim=c(range(dt$RT,na.rm=TRUE)), xlim=c(range(dt$TrialNo,na.rm=TRUE)),xaxs='i')
mod=rlm(dt$RT[hitpos]~dt$TrialNo[hitpos] )
abline(mod,lty=3,lwd=2, col='darkred')

lines(locpoly(dt$TrialNo[hitpos]  ,dt$RT[hitpos],bandwidth=RTbw), col='darkred')



plot(dt$DimmIntend[hitpos]  ,dt$ReleaseTime[hitpos], type='p',pch=20, col='red',ylim=c(range(dt$ReleaseTime,na.rm=TRUE)), xlim=c(range(dt$DimmIntend,na.rm=TRUE)),xaxs='i')
lines(locpoly(dt$DimmIntend[hitpos]  ,dt$ReleaseTime[hitpos],bandwidth=RTbw), col='darkred')


points(dt$DimmIntend[earlypos]  ,dt$ReleaseTime[earlypos], pch=20, col='green')
lines(locpoly(dt$DimmIntend[earlypos]  ,dt$ReleaseTime[earlypos],bandwidth=RTbw), col='darkgreen')

abline(a=0,b=1,,lty=3,lwd=2)



plot(dt$DimmOn[hitpos]  ,dt$RT[hitpos], type='p',pch=20, col='red',ylim=c(range(dt$RT,na.rm=TRUE)), xlim=c(range(dt$DimmIntend,na.rm=TRUE)),xaxs='i')
mod=rlm(dt$RT[hitpos]~dt$DimmOn[hitpos] )
abline(mod,lty=3,lwd=2, col='darkred')

lines(locpoly(dt$DimmOn[hitpos]  ,dt$RT[hitpos],bandwidth=RTbw), col='darkred')


