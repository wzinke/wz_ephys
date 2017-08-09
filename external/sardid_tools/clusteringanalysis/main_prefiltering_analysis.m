
clc;
clear;
close all;

suffix = 'prefilteringAnalysis';
if ~exist(suffix,'dir')
  mkdir(suffix)
end

% waveform Data
load('wFpreprocessed.mat');

% peak-to-trough duration, time of 25% repolarization and 1st PCA component of the two
var1 = wFpreprocessed.par(:,1);
var1 = var1(~isnan(wFpreprocessed.parPCA1stcomp));
var2 = wFpreprocessed.par(:,3);
var2 = var2(~isnan(wFpreprocessed.parPCA1stcomp));
var3 = wFpreprocessed.parPCA1stcomp(~isnan(wFpreprocessed.parPCA1stcomp));
dataWF = [var1,var2,var3];

% cell identification
datasetnameWF = wFpreprocessed.datasetname;
datasetnameWF = datasetnameWF(~isnan(wFpreprocessed.parPCA1stcomp));
spikechannelWF = wFpreprocessed.spikechannel;
spikechannelWF = spikechannelWF(~isnan(wFpreprocessed.parPCA1stcomp));
celltypeWF = wFpreprocessed.celltype(~isnan(wFpreprocessed.parPCA1stcomp));

dataWF_labels{1} = 'pk2trgh';
dataWF_labels{2} = 't4rplzn';
dataWF_labels{3} = '1CmpPCA';

% matrix of correlations
[dataWFcorr,dataWFcorr_pval]= corr(dataWF,'type','Spearman');

% clustering measures by correlation
Z = linkage(dataWF','single','spearman');
cutoff = 0.1; % cutoff: correlation of 0.9
figure('color','none','visible','off');
set(gca,'layer','top','color','none')
[h,t,outperm] = dendrogram(Z,'colorthreshold',cutoff,'labels',dataWF_labels);
% one way to decide on which measure to keep is the one that is closer to the center in a significant cluster
set(h,'color','k','linewidth',2);
xl = get(gca,'XTick');
hold on
plot([min(xl)-0.25 max(xl)+0.25],[cutoff cutoff],'r--','linewidth',2)
axis([min(xl)-0.25 max(xl)+0.25 0 1])
xlabel('Waveform measures','fontsize',16);
ylabel('Dissimilarity (Spearman correlation)','fontsize',16);
set(gca,'fontSize',16,'LineWidth',2,'TickDir','out','Box','off','YTick',0:.2:1)
plot2svg([suffix,'/wF.svg']);
close gcf;

datanonnormWF = dataWF(:,[1:2]);
% normalizalizing within [0,1]
datanormWF = standardization(dataWF(:,[1:2])); % 1stCompPCA is (by construction) redundant for clustering purposes

% activity data
Activity_Results = 'Activity.mat';
load(Activity_Results);

% measures
var1 = RES.FR;
var2 = RES.FF;
var3 = RES.Cv;
var4 = RES.Cv2;
var5 = RES.Lv;
var6 = RES.LvR;

dataACT = [var1(~isnan(var6)),var2(~isnan(var6)),var3(~isnan(var6)),var4(~isnan(var6)),var5(~isnan(var6)),var6(~isnan(var6))];

% cell identification
datasetnameACT = RES.datasetname;
datasetnameACT = datasetnameACT(~isnan(var6));
spikechannelACT = RES.spikechannel;
spikechannelACT = spikechannelACT(~isnan(var6));

for i = 1:length(datasetnameACT)
  tmp = datasetnameACT{i};
  tmp = tmp(1:findstr(tmp,'-')-1);
  datasetnameACT{i} = tmp;
end

clear RES;

dataACT_labels{1} = 'FR';
dataACT_labels{2} = 'FF';
dataACT_labels{3} = 'Cv';
dataACT_labels{4} = 'Cv2';
dataACT_labels{5} = 'Lv';
dataACT_labels{6} = 'LvR';

% matrix of correlations
[dataACTcorr,dataACTcorr_pval] = corr(dataACT,'type','Spearman');

% clustering measures by correlation
Z = linkage(dataACT','single','spearman');
cutoff = 0.1; % cutoff: correlation of 0.9
figure('color','none','visible','off');
set(gca,'layer','top','color','none')
[h,t,outperm] = dendrogram(Z,'colorthreshold',cutoff,'labels',dataACT_labels);
% one way to decide on which measure to keep is the one that is closer to the center in a significant cluster, thus, between Cv2, Lv and LvR, Lv is at center and then it should be kept while the other two should be discarded... unfortunately I didn't find a way to access this information without manually checking the dendrogram figure...
set(h,'linewidth',2);
xl = get(gca,'XTick');
hold on
plot([min(xl)-0.25 max(xl)+0.25],[cutoff cutoff],'r--','linewidth',2)
axis([min(xl)-0.25 max(xl)+0.25 0 1])
xlabel('Activity measures','fontsize',16);
ylabel('Dissimilarity (Spearman correlation)','fontsize',16);
set(gca,'fontSize',16,'LineWidth',2,'TickDir','out','Box','off','YTick',0:.2:1)
plot2svg([suffix,'/variability.svg']);
close gcf;

datanonnormACT = dataACT(:,[1:3,5]);
% normalizalizing within [0,1]
datanormACT = standardization(dataACT(:,[1:3,5]));

% putting the two datasets together for posterior clustering
cont = 0;
maxSize=max([size(datanormWF,1),size(datanormACT,1)]);
datanonnormAll=nan(maxSize,size(datanormWF,2)+size(datanormACT,2));
datanormAll=nan(maxSize,size(datanormWF,2)+size(datanormACT,2));
for i = 1:length(datasetnameWF)
  for j = 1:length(datasetnameACT)
    if strcmp(datasetnameWF{i},datasetnameACT{j}) && strcmp(spikechannelWF{i},spikechannelACT{j})
      cont = cont+1;
      datanonnormAll(cont,:)=[datanonnormWF(i,:),datanonnormACT(j,:)];
      datanormAll(cont,:)=[datanormWF(i,:),datanormACT(j,:)];
      celltype{cont}=celltypeWF{i};
      datasetname{cont}=datasetnameWF{i};
      spikechannel{cont}=spikechannelWF{i};
    end
  end
end
datanonnormAll = datanonnormAll(1:cont,:);
datanormAll = datanormAll(1:cont,:);
datanormAll_labels = {dataWF_labels{1:2},dataACT_labels{[1:3,5]}};

% variance contribution
var_datanormAll = var(datanormAll);
[sortedVar_datanormAll,indSortedVar_datanormAll] = sort(var_datanormAll,2,'descend');
explainedVar_datanormAll = 100*sortedVar_datanormAll/sum(sortedVar_datanormAll);
cumExplainedVar_datanormAll = cumsum(explainedVar_datanormAll);

% variance cutoff 90%, keeping only measurements that together account for >=90% of the variance
var_cutoff = 90;
indcutoff = find(cumExplainedVar_datanormAll >= var_cutoff,1);
datanonnorm4Cluster = datanonnormAll(:,indSortedVar_datanormAll);
datanonnorm4Cluster = datanonnorm4Cluster(:,1:indcutoff);
datanorm4Cluster = datanormAll(:,indSortedVar_datanormAll);
datanorm4Cluster = datanorm4Cluster(:,1:indcutoff);
datanorm4Cluster_labels = datanormAll_labels(indSortedVar_datanormAll);
datanorm4Cluster_labels = datanorm4Cluster_labels(1:indcutoff);
explainedVar_datanorm4Cluster = explainedVar_datanormAll(1:indcutoff);
cumExplainedVar_datanorm4Cluster = cumsum(explainedVar_datanorm4Cluster);

width = 1/3;
figure('color','none','visible','off');
hold on
set(gca,'layer','top','color','none')
h1=bar([1:6]-0.5*width,[explainedVar_datanorm4Cluster,100-cumExplainedVar_datanorm4Cluster(end)],'linewidth',2);
set(h1,'EdgeColor','k','FaceColor',[0.5 0.5 0.5],'BarWidth',width);
h2=bar([1:6]+0.5*width,[cumExplainedVar_datanorm4Cluster,100],'linewidth',2);
set(h2,'EdgeColor','k','FaceColor','w','BarWidth',width);
plot([0.5 6.5],[var_cutoff var_cutoff],'r--','linewidth',2)
xlabels = datanormAll_labels(indSortedVar_datanormAll);
legend([h1 h2],'Specific','Cumulative','fontSize',16,'Location',[.21,.65,.1,.1]); % left, bottom, width, height
legend boxoff;
ylabel('Variance explained (%)','fontSize',16)
axis([0.5 6.5 0 105])
set(gca,'yTick',0:20:100)
set(gca,'xTick',1:8)
set(gca,'XTickLabel',xlabels)
set(gca,'fontSize',16,'LineWidth',2,'TickDir','out','Box','off')
filename=[suffix,'/variance_explained.svg'];
plot2svg(filename);

[datanormClusterCorr,datanormClusterCorr_pval] = corr(datanorm4Cluster,'type','Spearman');

% save data
save([suffix,'/clusterDataset.mat'],'datasetname','spikechannel','datanonnorm4Cluster','datanorm4Cluster','datanorm4Cluster_labels','explainedVar_datanorm4Cluster','cumExplainedVar_datanorm4Cluster','datanormClusterCorr','datanormClusterCorr_pval','celltype','-v7');
