
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

% control cluster data, checked clusters from 5:15 (offset = 4), the relevant number of clusters is 7
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

% run randomizations
randomizedClusters.labels = nan(numRand,nElements);
randomizedClusters.randPermLabels = nan(numRand,nElements);
% random selection of neurons, keeping the sample size and allowing for repetitions
randElems = sort(randi(nElements,numRand,nElements),2);
% centers
randomizedClusters.centers = nan(numRand,relNumCluster,numDim);
randomizedClusters.randPermCenters = nan(numRand,relNumCluster,numDim);
% distances
randomizedClusters.distances = nan(numRand,relNumCluster,relNumCluster);
randomizedClusters.randPermDistances = nan(numRand,relNumCluster,relNumCluster);
% cluster to control assignment
randomizedClusters.cluster2ControlAssoc = nan(numRand,relNumCluster);
randomizedClusters.randPermCluster2ControlAssoc = nan(numRand,relNumCluster);
% percentage of matched elements in clusters with respect to control
randomizedClusters.percentMatch = nan(numRand,relNumCluster);
randomizedClusters.randPermPercentMatch = nan(numRand,relNumCluster);

for i = 1:numRand
  % obtaining the metaclusters... :
  display('running metacluster:');
  display(i);
  tmpData = control.data(randElems(i,:),:);
  randomizedClusters.labels(i,:) = function_pairingOfClusterElements(type,tmpData);

  % random permutation of cluster elements through labels
  randomizedClusters.randPermLabels(i,:) = control.labels(randperm(nElements));

  % compute centers
  for j = 1:relNumCluster
    randomizedClusters.centers(i,j,:) = mean(tmpData(randomizedClusters.labels(i,:)==j,:),1);
    randomizedClusters.randPermCenters(i,j,:) = mean(control.data(randomizedClusters.randPermLabels(i,:)==j,:),1);

    % compute distances with respect to control clusters
    for k = 1:relNumCluster
      randomizedClusters.distances(i,j,k) = pdist2(squeeze(randomizedClusters.centers(i,j,:))',control.centers(k,:));
      randomizedClusters.randPermDistances(i,j,k) = pdist2(squeeze(randomizedClusters.randPermCenters(i,j,:))',control.centers(k,:));
    end
  end

  % assign clusters: closest to each control cluster
  for j = 1:relNumCluster
    [~,randomizedClusters.cluster2ControlAssoc(i,j)] = min(randomizedClusters.distances(i,:,j));
    [~,randomizedClusters.randPermCluster2ControlAssoc(i,j)] = min(randomizedClusters.randPermDistances(i,:,j));

    % compute percentage of unique elements that match with respect to control clusters
    randElemInCluster = randElems(i,randomizedClusters.labels(i,:)==randomizedClusters.cluster2ControlAssoc(i,j));
    randomizedClusters.percentMatch(i,j) = 100*length(intersect(randElemInCluster,control.elemInCluster{j}))/length(intersect(randElems(i,:),control.elemInCluster{j}));
    randPermElemInCluster = find(randomizedClusters.randPermLabels(i,:)==randomizedClusters.randPermCluster2ControlAssoc(i,j));
    randomizedClusters.randPermPercentMatch(i,j) = 100*length(intersect(randPermElemInCluster,control.elemInCluster{j}))/length(control.elemInCluster{j});
  end
end

% for each control cluster compute:
% 1. mean and ste of the percentage of match
randomizedClusters.meanPercentMatch = mean(randomizedClusters.percentMatch);
randomizedClusters.stePercentMatch = ste(randomizedClusters.percentMatch);
randomizedClusters.meanRandPermPercentMatch = mean(randomizedClusters.randPermPercentMatch);
randomizedClusters.steRandPermPercentMatch = ste(randomizedClusters.randPermPercentMatch);

width = 1/3;
figure('color','none','visible','off');
hold on
set(gca,'layer','top','color','none')
h1=bar((1:relNumCluster)-.5*width,randomizedClusters.meanPercentMatch,'linewidth',1);
set(h1,'EdgeColor','k','FaceColor',[0.6 0.6 0.6],'BarWidth',width);
plot(ones(2,1)*(1:relNumCluster)-.5*width,[randomizedClusters.meanPercentMatch+randomizedClusters.stePercentMatch;randomizedClusters.meanPercentMatch-randomizedClusters.stePercentMatch],'k','linewidth',2);
h2=bar((1:relNumCluster)+.5*width,randomizedClusters.meanRandPermPercentMatch,'linewidth',1);
set(h2,'EdgeColor','k','FaceColor','w','BarWidth',width);
plot(ones(2,1)*(1:relNumCluster)+.5*width,[randomizedClusters.meanRandPermPercentMatch+randomizedClusters.steRandPermPercentMatch;randomizedClusters.meanRandPermPercentMatch-randomizedClusters.steRandPermPercentMatch],'k','linewidth',2);
plot([0,relNumCluster+1],100/relNumCluster*[1 1],'r--','linewidth',1)
xlabel('Cell classes','fontSize',16)
ylabel('% of matching cells','fontSize',16)
axis([0.5 relNumCluster+0.5 0 100])
legend([h1 h2],'Randomized clustering','Random assignments','fontSize',16,'Location','NorthOutside');
legend boxoff;
set(gca,'xTick',1:relNumCluster,'yTick',0:20:100)
set(gca,'XTickLabel',{'B1','B2','B3','B4','N1','N2','N4'})
set(gca,'fontSize',16,'LineWidth',1,'TickDir','out','Box','off')
filename=[figsdir,'/percMatchedCells.svg'];
plot2svg(filename);

% 2. difference of not assoc minus assoc cluster distances: mean and ste
randomizedClusters.distAssocClusters = nan(numRand,relNumCluster);
randomizedClusters.randPermDistAssocClusters = nan(numRand,relNumCluster);
randomizedClusters.distNotAssocClusters = nan(numRand,relNumCluster);
randomizedClusters.randPermDistNotAssocClusters = nan(numRand,relNumCluster);
for i = 1:numRand
  tmp1 = nan(relNumCluster,relNumCluster-1);
  tmp2 = nan(relNumCluster,relNumCluster-1);
  for j = 1:relNumCluster
    randomizedClusters.distAssocClusters(i,j) = randomizedClusters.distances(i,randomizedClusters.cluster2ControlAssoc(i,j),j);
    randomizedClusters.randPermDistAssocClusters(i,j) = randomizedClusters.randPermDistances(i,randomizedClusters.randPermCluster2ControlAssoc(i,j),j);
    tmp1(j,:) = randomizedClusters.distances(i,setdiff(1:relNumCluster,randomizedClusters.cluster2ControlAssoc(i,j)),j);
    tmp2(j,:) = randomizedClusters.randPermDistances(i,setdiff(1:relNumCluster,randomizedClusters.randPermCluster2ControlAssoc(i,j)),j);
  end
  randomizedClusters.distNotAssocClusters(i,:) = nanmean(tmp1,2);
  randomizedClusters.randPermDistNotAssocClusters(i,:) = mean(tmp2,2);
end
randomizedClusters.meanDiffDistance = mean(randomizedClusters.distNotAssocClusters-randomizedClusters.distAssocClusters,1);
randomizedClusters.steDiffDistance = ste(randomizedClusters.distNotAssocClusters-randomizedClusters.distAssocClusters);
randomizedClusters.meanDiffRandPermDistance = mean(randomizedClusters.randPermDistNotAssocClusters-randomizedClusters.randPermDistAssocClusters,1);
randomizedClusters.steDiffRandPermDistance = ste(randomizedClusters.randPermDistNotAssocClusters-randomizedClusters.randPermDistAssocClusters);

width = 1/3;
figure('color','none','visible','off');
hold on
set(gca,'layer','top','color','none')
h1=bar((1:relNumCluster)-.5*width,randomizedClusters.meanDiffDistance,'linewidth',1);
set(h1,'EdgeColor','k','FaceColor',[0.6 0.6 0.6],'BarWidth',width);
plot(ones(2,1)*(1:relNumCluster)-.5*width,[randomizedClusters.meanDiffDistance+randomizedClusters.steDiffDistance;randomizedClusters.meanDiffDistance-randomizedClusters.steDiffDistance],'k','linewidth',2);
h2=bar((1:relNumCluster)+.5*width,randomizedClusters.meanDiffRandPermDistance,'linewidth',1);
set(h2,'EdgeColor','k','FaceColor','w','BarWidth',width);
plot(ones(2,1)*(1:relNumCluster)+.5*width,[randomizedClusters.meanDiffRandPermDistance+randomizedClusters.steDiffRandPermDistance;randomizedClusters.meanDiffRandPermDistance-randomizedClusters.steDiffRandPermDistance],'k','linewidth',2);
xlabel('Cell classes','fontSize',16)
ylabel('Cluster distance (not assoc. vs assoc.)','fontSize',16)
xlim([0.5 relNumCluster+0.5])
legend([h1 h2],'Randomized clustering','Random assignments','fontSize',16,'Location','NorthOutside');
legend boxoff;
set(gca,'xTick',1:relNumCluster,'yTick',0:0.1:0.5)
set(gca,'XTickLabel',{'B1','B2','B3','B4','N1','N2','N4'})
set(gca,'fontSize',16,'LineWidth',1,'TickDir','out','Box','off')
filename=[figsdir,'/diffDistanceBetweenClusters.svg'];
plot2svg(filename);

% total elapsed time in the analysis
elapsedTime=toc(timeId);

elapsedTimeHours=floor(elapsedTime/3600);
elapsedTimeMinuts=floor(mod(elapsedTime,3600)/60);
elapsedTimeSeconds=mod(mod(elapsedTime,3600),60);

fprintf('\n>>>> total time = %d hours, %d minuts, %.3f seconds\n\n', elapsedTimeHours,elapsedTimeMinuts,elapsedTimeSeconds);
