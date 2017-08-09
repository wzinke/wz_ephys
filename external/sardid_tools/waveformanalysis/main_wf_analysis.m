clear;
close all;
clc;

FIGSDIR = 'waveformAnalysisPanels';
if ~exist(FIGSDIR,'dir')
    mkdir(FIGSDIR);
end

%%%% loading waveform data %%%%
load wForig.mat;

%%%% preprocessing waveforms %%%%

% setting up configuration for waveForm preprocessing: trough alignment, interpolation, normalization, and sampling frequency
cfg = struct('troughalign',1,'interpn',10,'interptype','splines','normalize',1,'fsample',40000);

[W, X, par, parName] = waveformPreprocessing(wForig.W,cfg,wForig.isolationquality); % isolationquality: 1=multi-unit, 2=mostly-single-unit, 3=single-unit, 0=not-specified

wFpreprocessed.cfg = cfg;
wFpreprocessed.W = W;
wFpreprocessed.X = X;
wFpreprocessed.par = par;
wFpreprocessed.parName = parName;
wFpreprocessed.datasetname = wForig.datasetname;
wFpreprocessed.spikechannel = wForig.spikechannel;
wFpreprocessed.Worig = wForig.W;
wFpreprocessed.Xorig = wForig.X;
wFpreprocessed.isolationquality = wForig.isolationquality;
wFpreprocessed.isolationqualityInfo = wForig.isolationqualityInfo;

%%%% separating cell types from waveform measurements %%%%

% index to highly isolated units
sel = find(wFpreprocessed.isolationquality==3);

%%%% par #1 peak-to-trough duration %%%%
X1 = par(sel,1);
X1 = X1(~isnan(X1));

mu = [.17;.35]; % initial mean of each component in the 2-Gaussian fit
sigma = [.02^2;.02^2]; % initial covariance matrix of each component in the 2-Gaussian fit
pcomponents = [.1;.9]; % initial mixing proportions of each component in the 2-Gaussian fit
Xlim = [0.1 0.68]; % xlim for plotting

[dip1,pdip1,xl1,xu1,ind_nar1,ind_bro1,ind_fuz1,aic1_1,aic1_2,bic1_1,bic1_2] = waveformSeparation(X1,'peak2trough','Peak to trough duration',Xlim,mu,sigma,pcomponents,FIGSDIR);

%%%% par #3 duration of (25%) repolarization %%%%
X3 = par(sel,3);
X3 = X3(~isnan(X3));

mu = [.13;.28];
sigma = [.05^2;.1^2];
pcomponents = [.3;.7];
Xlim = [0.03 0.5];
[dip3,pdip3,xl3,xu3,ind_nar3,ind_bro3,ind_fuz3,aic3_1,aic3_2,bic3_1,bic3_2] = waveformSeparation(X3,'25percRepol','Time for repolarization',Xlim,mu,sigma,pcomponents,FIGSDIR);

%%%% histogram for the first component of the PCA for par #1 and par #3 %%%%

Xpc1 = par(sel,1);
Xpc2 = par(sel,3);
Ipc = find(~isnan(Xpc2));
Xpc1 = Xpc1(Ipc);
Xpc2 = Xpc2(Ipc);

% how correlated are the two measures
[rcorr,pcorr] = corr(Xpc1,Xpc2,'type','Pearson');

% PCA
Xpca = [Xpc1 Xpc2];
[coeff,score] = pca(Xpca,'Economy','on');

% data reconstruction (takes into account whether data is centered):
% Xpca = score*coeff'+ones(size(Xpca-ones,1),1)*mean(Xpca,1);
% as inv(coeff) = coeff', then
% score = (Xpca-ones(size(Xpca-ones,1),1)*mean(Xpca,1))*coeff;

latent = var(score);
explained = 100*latent/sum(latent);

% first component of the PCA
Xpca1stcomp = score(:,1);

numUnits = size(W,1);
wFpreprocessed.parPCA1stcomp = nan(numUnits,1);
wFpreprocessed.parPCA1stcomp(sel(Ipc)) = Xpca1stcomp;

% variance explained by the first component of the PCA
FirstPCAvarExplained = explained(1);

mu = [-.17;.07];
sigma = [.04^2;.06^2];
pcomponents = [.2;.8];
Xlim = [-0.3 0.35];
[dipPCA,pdipPCA,xlPCA,xuPCA,ind_narPCA,ind_broPCA,ind_fuzPCA,aicPCA_1,aicPCA_2,bicPCA_1,bicPCA_2] = waveformSeparation(Xpca1stcomp,'PCA1stComp','First component of the PCA',Xlim,mu,sigma,pcomponents,FIGSDIR);

% plotting waveforms

toi = find(X>=-0.5 & X<=1.0);

figure('color','none','visible','off');
hold on
set(gca,'layer','top','color','none')
plot([min(X(toi)) max(X(toi))],[1 1],'k-','linewidth',1)
plot(X(toi),W(sel(Ipc(ind_broPCA)),toi)','b-','linewidth',1)
plot(X(toi),W(sel(Ipc(ind_fuzPCA)),toi)','color',[0.1 0.1 0.1],'linewidth',1)
plot(X(toi),W(sel(Ipc(ind_narPCA)),toi)','r-','linewidth',1)
axis([-0.5,1,-1,1.05])
xlabel('Time (ms)','fontSize',16)
ylabel('Normalized voltage','fontSize',16)
set(gca,'fontSize',16,'LineWidth',1,'TickDir','out','Box','off','YTick',-1:0.5:1)
plot2svg([FIGSDIR,'/waveforms.svg'])

% example waveform prototype:
ind_prototype=77;
% open('prototypeWF.fig');
% y = get(gco,'ydata');
% ind_prototype = strmatch(num2str(y),num2str(W(sel(Ipc(ind_nar)),toi)));

figure('color','none','visible','off');
hold on
set(gca,'layer','top','color','none')
plot(X(toi),W(sel(Ipc(ind_narPCA(ind_prototype))),toi),'k-','linewidth',2)
ylim([-1.05,1.05])
axis off
plot2svg([FIGSDIR,'/prototypeWF.svg'])

figure('color','none','visible','off');
hold on
set(gca,'layer','top','color','none')
plot(X(toi),W(sel(Ipc(ind_narPCA(ind_prototype))),toi),'ko','linewidth',2)
ylim([-1.05,1.05])
axis off
plot2svg([FIGSDIR,'/prototypeWF_interpolated.svg'])

figure('color','none','visible','off');
hold on
set(gca,'layer','top','color','none')
plot(wForig.X,wForig.W(sel(Ipc(ind_narPCA(ind_prototype))),:),'ko','linewidth',2)
xlim([-.05,1.325])
axis off
plot2svg([FIGSDIR,'/prototypeWF_notInterpolated.svg'])

% tagging celltypes

wFpreprocessed.celltype = cell(1,length(wFpreprocessed.isolationquality));
wFpreprocessed.celltype(:) = {'not-defined'};
wFpreprocessed.celltype(sel(Ipc(ind_narPCA))) = {'narrow'};
wFpreprocessed.celltype(sel(Ipc(ind_broPCA))) = {'broad'};
wFpreprocessed.celltype(sel(Ipc(ind_fuzPCA))) = {'fuzzy'};

% New: post-hoc analysis for those units for which we do not have the time for repolarization; we will consider 'broad' spiking units those with peak-to-trough durations longer than percentile 5 of all BS, and 'narrow' those shorter than percentile 95 of all NS, the rest will remain unclassified: 'fuzzy'

cutOffBS = quantile(Xpc1(ind_broPCA),.05);
cutOffNS = quantile(Xpc1(ind_narPCA),.95);

wFpreprocessed.celltype(sel(isnan(par(sel,3)) & par(sel,1) >= cutOffBS)) = {'broad'};
wFpreprocessed.celltype(sel(isnan(par(sel,3)) & par(sel,1) <= cutOffNS)) = {'narrow'};
wFpreprocessed.celltype(sel(isnan(par(sel,3)) & par(sel,1) > cutOffNS & par(sel,1) < cutOffBS)) = {'fuzzy'};

save('wFpreprocessed.mat','wFpreprocessed','-v7');
