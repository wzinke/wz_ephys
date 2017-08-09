
clear;
close all;
clc;

% to control the time of the analysis
timeId=tic;

numRand = 200;
figsdir = 'validationAnalysis';

% load control data

type = 'Kmeans';
resultsdir = 'clusteringResults';
matfile = [resultsdir,'/',type,'ClusteringResults'];

% control cluster data checked clusters from 5:10 (offset = 4), the relevant number of clusters is 7
offset = 4;
relNumCluster = 7;
relIndex = relNumCluster-offset;
tmp = load(matfile);
control.data = tmp.datanorm4Cluster;
numDim = size(control.data,2);
nElements = length(tmp.labels);
control.labels = nan(1,nElements);
control.centers = nan(relNumCluster,numDim);
control.clustOrder = tmp.clustOrder;
cont = 0;
for i = control.clustOrder
  if (i == control.clustOrder(7)) % filtering out the non-reliable cluster
    continue
  else
    cont = cont + 1; % cont 1:7, representing [B1:B4,N1,N2,N4], N3 is not reliable
    control.labels(tmp.labels(relIndex,:)==i) = cont;
    control.centers(cont,:) = mean(control.data(control.labels==cont,:),1);
    control.elemInCluster{cont} = tmp.clustFilt{relIndex,i};
  end
end
clear tmp;

% monkeys info
infodir = 'prefilteringAnalysis';
matfile = [infodir,'/clusterDataset'];
load(matfile,'datasetname');
mi.elems = find(strncmp(datasetname,'mi',2));
ry.elems = find(strncmp(datasetname,'ry',2));
% obtaining the metaclusters... :
display('running metacluster for monkeys:');
mi.data = control.data(mi.elems,:);
mi.labels = function_pairingOfClusterElements(type,mi.data);
ry.data = control.data(ry.elems,:);
ry.labels = function_pairingOfClusterElements(type,ry.data);
% centers
mi.centers = nan(relNumCluster,numDim);
ry.centers = nan(relNumCluster,numDim);
% distances
mi.distances = nan(relNumCluster,relNumCluster);
ry.distances = nan(relNumCluster,relNumCluster);
% cluster to control assignment
mi.cluster2ControlAssoc = nan(1,relNumCluster);
ry.cluster2ControlAssoc = nan(1,relNumCluster);
% percentage of matched elements in clusters with respect to control
mi.percentMatch = nan(1,relNumCluster);
ry.percentMatch = nan(1,relNumCluster);
for j = 1:relNumCluster
  mi.centers(j,:) = mean(mi.data(mi.labels==j,:),1);
  ry.centers(j,:) = mean(ry.data(ry.labels==j,:),1);
  % compute distances with respect to control clusters
  for k = 1:relNumCluster
    mi.distances(j,k) = pdist2(squeeze(mi.centers(j,:)),control.centers(k,:));
    ry.distances(j,k) = pdist2(squeeze(ry.centers(j,:)),control.centers(k,:));
  end
end
% assign clusters: closest to each control cluster and the percentage of match
for j = 1:relNumCluster
  [~,mi.cluster2ControlAssoc(j)] = min(mi.distances(:,j));
  [~,ry.cluster2ControlAssoc(j)] = min(ry.distances(:,j));
  % compute percentage of unique elements that match with respect to control clusters
  elemInCluster = mi.elems(mi.labels==mi.cluster2ControlAssoc(j));
  mi.percentMatch(j) = 100*length(intersect(elemInCluster,control.elemInCluster{j}))/length(intersect(mi.elems,control.elemInCluster{j}));
  elemInCluster = ry.elems(ry.labels==ry.cluster2ControlAssoc(j));
  ry.percentMatch(j) = 100*length(intersect(elemInCluster,control.elemInCluster{j}))/length(intersect(ry.elems,control.elemInCluster{j}));
end

width = 1/3;
figure('color','none','visible','off');
hold on
set(gca,'layer','top','color','none')
h1=bar((1:relNumCluster)-.5*width,mi.percentMatch,'linewidth',1);
set(h1,'EdgeColor','k','FaceColor','r','BarWidth',width);
h2=bar((1:relNumCluster)+.5*width,ry.percentMatch,'linewidth',1);
set(h2,'EdgeColor','k','FaceColor','b','BarWidth',width);
plot([0,relNumCluster+1],100/relNumCluster*[1 1],'k--','linewidth',1)
xlabel('Cell classes','fontSize',16)
ylabel('% of matching cells','fontSize',16)
axis([0.5 relNumCluster+0.5 0 100])
legend([h1 h2],'Monkey M','Monkey R','fontSize',16,'Location','NorthOutside');
legend boxoff;
set(gca,'xTick',1:relNumCluster,'yTick',0:20:100)
set(gca,'XTickLabel',{'B1','B2','B3','B4','N1','N2','N4'})
set(gca,'fontSize',16,'LineWidth',1,'TickDir','out','Box','off')
filename=[figsdir,'/percMatchedCells_monkeys.svg'];
plot2svg(filename);

% 2. difference of not assoc minus assoc cluster distances: mean and ste
mi.distAssocClusters = nan(1,relNumCluster);
ry.distAssocClusters = nan(1,relNumCluster);
mi.distNotAssocClusters = nan(1,relNumCluster);
ry.distNotAssocClusters = nan(1,relNumCluster);

tmp1 = nan(relNumCluster,relNumCluster-1);
tmp2 = nan(relNumCluster,relNumCluster-1);

for j = 1:relNumCluster
  mi.distAssocClusters(j) = mi.distances(mi.cluster2ControlAssoc(j),j);
  tmp1(j,:) = mi.distances(setdiff(1:relNumCluster,mi.cluster2ControlAssoc(j)),j);
  ry.distAssocClusters(j) = ry.distances(ry.cluster2ControlAssoc(j),j);
  tmp2(j,:) = ry.distances(setdiff(1:relNumCluster,ry.cluster2ControlAssoc(j)),j);
end
mi.distNotAssocClusters = nanmean(tmp1,2)';
ry.distNotAssocClusters = nanmean(tmp2,2)';
mi.diffDistance = mi.distNotAssocClusters-mi.distAssocClusters;
ry.diffDistance = ry.distNotAssocClusters-ry.distAssocClusters;

width = 1/3;
figure('color','none','visible','off');
hold on
set(gca,'layer','top','color','none')
h1=bar((1:relNumCluster)-.5*width,mi.diffDistance,'linewidth',1);
set(h1,'EdgeColor','k','FaceColor','r','BarWidth',width);
h2=bar((1:relNumCluster)+.5*width,ry.diffDistance,'linewidth',1);
set(h2,'EdgeColor','k','FaceColor','b','BarWidth',width);
xlabel('Cell classes','fontSize',16)
ylabel('Cluster distance (not assoc. vs assoc.)','fontSize',16)
xlim([0.5 relNumCluster+0.5])
legend([h1 h2],'Monkey M','Monkey R','fontSize',16,'Location','NorthOutside');
legend boxoff;
set(gca,'xTick',1:relNumCluster,'yTick',0:0.1:0.5)
set(gca,'XTickLabel',{'B1','B2','B3','B4','N1','N2','N4'})
set(gca,'fontSize',16,'LineWidth',1,'TickDir','out','Box','off')
filename=[figsdir,'/diffDistanceBetweenClusters_monkeys.svg'];
plot2svg(filename);

% total elapsed time in the analysis
elapsedTime=toc(timeId);

elapsedTimeHours=floor(elapsedTime/3600);
elapsedTimeMinuts=floor(mod(elapsedTime,3600)/60);
elapsedTimeSeconds=mod(mod(elapsedTime,3600),60);

fprintf('\n>>>> total time = %d hours, %d minuts, %.3f seconds\n\n', elapsedTimeHours,elapsedTimeMinuts,elapsedTimeSeconds);
