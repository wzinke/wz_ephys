
%%%%%%%%%%%%%%%%%%%
% Figure 1 Panels %
%%%%%%%%%%%%%%%%%%%


load([pwd,'/Workspace/GlauberIsingData.mat']);
tauIndex=2;
statename = '\tau=2';
BinIndex=1;
BinSize=10;
TempTab = [1 2 3];
TempList = (1:1:5);
tauGraph = [3 5];
TempGraph = [2 3];

nbV=2;
nbH=2;
spaceH=0.18;
spaceV=0.1;
borderH=0.06;
borderV=0.06;
HGraph = (1-(nbH-1)*spaceH-2*borderH)/nbH;
VGraph = (1-(nbV-1)*spaceV - 2*borderV)/nbV;

FigsParams;


Fig9Preds(PexpTab,PthTab,PthNoJ1Tab,PthNoJTab,FeaturesTab,BinIndex,BinSize,TempTab,tauIndex,datalen)
