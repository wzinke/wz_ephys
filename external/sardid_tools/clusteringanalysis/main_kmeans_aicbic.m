
clear;
close all;
clc;

% see http://stats.stackexchange.com/questions/55147/k-means-bic-to-validate-clusters-in-r

suffix = 'clusteringResults';
type = 'Kmeans';
% control cluster data checked clusters from 5:10 (offset = 4), the relevant number of clusters is 7
offset = 4;
minNumOfElementsCluster = 5;
for k = 5:15 % num of clusters
  relIndex = k-offset;
  tmp = load([suffix,'/',type,'ClusteringResults']);
  control.data = tmp.datanorm4Cluster;
  numDim = size(control.data,2);
  nElements = length(tmp.labels);
  control.labels = nan(1,nElements);
  control.clustOrder = tmp.clustOrder;
  cont = 0;
  for i = control.clustOrder
    if (length(tmp.clustFilt{relIndex,i}) < minNumOfElementsCluster)
     continue
    else
      cont = cont + 1; % cont 1:7, representing [B1:B4,N1,N2,N4], N3 is not reliable
      control.labels(tmp.labels(relIndex,:)==i) = cont;
      control.centers(cont,:) = mean(control.data(control.labels==cont,:),1);
      control.elemInCluster{cont} = tmp.clustFilt{relIndex,i};
      Nc(cont) = length(control.elemInCluster{cont});
      Vc(:,cont) = var(control.data(control.labels==cont,:),[],1); % normalizing by n-1, there is no cluster with only 1 element
    end
  end
  clear tmp;

  numRelClusters = cont;

  % V = var(control.data,[],1)'*ones(1,numRelClusters); % normalizing by n-1, there is no cluster with only 1 element
  % logL = -Nc.*sum(.5*log(Vc+V),1);
  logL = -Nc.*sum(.5*log(Vc),1); % there is no cluster with only 1 element
  aic(relIndex) = -2 * sum(logL) + 4*numRelClusters*numDim;
  bic(relIndex) = -2 * sum(logL) + 2*numRelClusters*numDim*log(nElements);
end
display(aic);
display(bic);

