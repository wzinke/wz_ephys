
clear;
close all;
clc;

type = 'Kmeans';
resultsdir = 'clusteringResults';
matfile = [resultsdir,'/',type,'ClusteringResults'];
load(matfile);

% calculating centers of k-means clusters
numClusters = 7;
offset = 4; % it starts with numClusters = 5
numOfDim=size(datanorm4Cluster,2);

nC=length(clustFilt(1,:));
k=0;
tmp=1:size(datanorm4Cluster,1);
for i=1:nC
    if(~isempty([clustFilt{numClusters-offset,i}]))
        k = k+1;
        centers(k,:) = mean(datanorm4Cluster([clustFilt{numClusters-offset,i}],:),1);
        centersSD(k,:) = std(datanorm4Cluster([clustFilt{numClusters-offset,i}],:),0,1);
        numElem(k) = length([clustFilt{numClusters-offset,i}]);
        tmp = setdiff(tmp,[clustFilt{numClusters-offset,i}]);
    end
end

% remaining elements if any
if ~isempty(tmp)
  k = k+1;
  clustFilt{numClusters-offset,k} = tmp;
  numElem(k) = size(datanorm4Cluster,1)-sum(numElem);
  centers(k,:) = mean(datanorm4Cluster([clustFilt{numClusters-offset,k}],:),1);
  centersSD(k,:) = std(datanorm4Cluster([clustFilt{numClusters-offset,k}],:),0,1);
end

Ap = [];
for i = 1:k
  % I add a term only to keep the proper structure of the dendrogram according to the neurons in the clusters
  Ap = [Ap;(ones(numElem(i),1)+0.01*((1:numElem(i))'/numElem(i)-0.5))*centers(i,:)];
end
Y = pdist(Ap);
Z = linkage(Y,'average');
c = cluster(Z,'maxclust',k);

% threshold color
t = sort(Z(:,3));
Cth = t(size(Z,1)+2-k);

figure('color','none','visible','off');
set(gca,'layer','top','color','none')
[h,~,outperm] = dendrogram(Z,0,'ColorThreshold',Cth,'orientation','left');
hold on
set(h,'LineWidth',1)
axis off
ylabel('Neurons','fontsize',16);
xlabel('Distance','fontsize',16);
set(gca,'fontSize',16,'LineWidth',1,'TickDir','out','Box','off','XTick',0:.2:1,'YTick',[])
plot2svg([resultsdir,'/dendrogram_',type,'_centroids.svg'])

% matching the dendrogram with the results of the kmeans to determine the cluster ordering
clustersNumOfElements = find(diff(c(outperm))~=0);
clustersNumOfElements = [clustersNumOfElements(1);diff([clustersNumOfElements;length(c)])];
for i = 1:length(clustersNumOfElements)
  clustOrder(i) = find(clustersNumOfElements(i)==numElem);
end

save(matfile,'clustFilt','clustOrder','-append');
