# setwd('/DATA/Ephys/TST_data/tst')
#
# flnm = '/DATA/Ephys/TST_data/src/Q/TL/Q062305002_DSP04a_R_FEF_SEARCH_SPIKE.mat'
# timwin = c(40,160)

#https://stat.ethz.ch/pipermail/r-help/2006-March/102204.html

WZ_getRespDiff = function(flnm, timwin = c(40,160)

require('R.matlab')
require('R.matlab')

dt = readMat(flnm)

p = dt$correct ==1 & dt$SetSZ == 4 & dt$Tpos != 255 & is.na(dt$MStim.[,1])

spkrts = dt$spikedensities[p,]
Tpos   = as.factor(dt$Tpos[p])

ntm = dim(spkrts)[2]
pvec = rep(NA,ntm)

trialno = seq(dt$nTrials)[p]

for(i in seq(ntm){

    ctm = spkrts[,i]
    aovmod = lmRob(ctm~trialno+Tpos, na.action=na.exclude)

    aovmod = aov(ctm~trialno+Tpos)
    kruskal.test(ctm~trialno)
}