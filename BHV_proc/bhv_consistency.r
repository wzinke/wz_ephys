bhvdir = '/SCRATCH/BHV_files'
tblfl  = 'Ingrid_get_joystick_2014_04_26.dat'
bhvfl  = 'Ingrid_2014_37_26_trialtable.dat'
# bhvfl  = 'Ingrid_2014_59_27_trialtable.dat'

setwd(bhvdir)

tbl = read.table(tblfl, header=TRUE)
bhv = read.table(bhvfl, header=TRUE, sep=',')

tbl$Tstart  = tbl$Tstart  * 1000;
tbl$TFixOn  = tbl$TFixOn  * 1000 - tbl$Tstart 
tbl$DimmEff = tbl$DimmEff * 1000 - tbl$Tstart 
tbl$TRel    = tbl$TRel    * 1000 - tbl$Tstart 

plot(abs(tbl$TFixOn),tbl$WaitStart)

library(MASS)

plot(abs(tbl$TFixOn),abs(tbl$TFixOn)-tbl$WaitStart,xlab='wait start [ms]',ylab='time difference [ms]',pch=20)
mod=rlm((abs(tbl$TFixOn)-tbl$WaitStart)~bhv$WaitStart)
abline(mod,lty=3,lwd=2,col='darkred')


truehist(abs(tbl$TFixOn)-tbl$WaitStart,col='black',h=1,prob=FALSE,xlab='time difference [ms]',ylab='count',main='toc vs. trialtime')


truehist(bhv$WaitStart -tbl$WaitStart,col='black',h=1,prob=FALSE,xlab='time difference [ms]',ylab='count',main='event vs. trialtime')



truehist(bhv$WaitStart -abs(tbl$TFixOn),col='black',h=1,prob=FALSE,xlab='time difference [ms]',ylab='count',main='event vs. toc')

# reward lag
tblRew=tbl$RewTime-tbl$ReleaseTime
truehist(tblRew,col='black',border='grey',h=1,prob=FALSE,xlab='time [ms]',ylab='count',main='reward delay')

truehist(bhv$WaitRew,col='black',border='grey',h=1,prob=FALSE,xlab='time [ms]',ylab='count',main='reward delay')

truehist(bhv$RewLagEff,col='black',border='grey',h=1,prob=FALSE,xlab='time [ms]',ylab='count',main='reward delay')
truehist(bhv$RewLag,col='black',border='grey',breaks=100.05+c(-2,-1,0,1,2,3,4,5,6,7,8)/5,prob=FALSE,xlab='time [ms]',ylab='count',main='reward delay')


truehist(tblRew-bhv$WaitRew,col='black',border='grey',prob=FALSE,xlab='time [ms]',ylab='count',main='Reward Consistency')


# RT
hpos = bhv$TrialError == 0
RTbw=50
truehist(bhv$RT[hpos],col='black',border='grey',h=RTbw,prob=TRUE,xlab='time [ms]',ylab='probability',main='RT: eyejoytrack')
lines(density(bhv$RT[hpos],  bw=RTbw, from=0, cut=0, na.rm=TRUE),col='green',lwd=2.5)

truehist(bhv$RT[hpos]-bhv$RTcalc[hpos],col='black',border='grey',h=0.25,prob=TRUE,xlab='time [ms]',ylab='probability',main='RT difference: eyejoytrack - events')


RTtbl = (tbl$TRel - tbl$DimmEff)
truehist(tbl$RT[hpos]-RTtbl[hpos],col='black',border='grey',h=0.25,prob=TRUE,xlab='time [ms]',ylab='probability',main='RT difference: events - toc')

truehist(bhv$RTcalc[hpos]-RTtbl[hpos],col='black',border='grey',h=0.25,prob=TRUE,xlab='time [ms]',ylab='probability',main='RT difference: eyejoytrack - toc')
