
clear;
close all;
clc;

% This script runs different realizations of the Kmeans/Kmedians clustering
% algorithm for different number of clusters to evaluate
% what number could be best

type = 'Kmeans';
% type = 'Kmedians';
if strcmp(type,'Kmeans')
  distance = 'sqEuclidean';
elseif strcmp(type,'Kmedians')
  distance = 'cityblock';
else
  error('type of clustering has to be either Kmeans or Kmedians');
end

% to control the time of the analysis
timeId=tic;

% number of clusters to consider
minNumberOfClusters=1;
maxNumberOfClusters=40;

numberOfRepeats=10; % number of realizations (to compare behavior for different random seeds)
numOfReplicates=50; % number of replicates (to then get the one with lowest energy)

nk=1; % skip computation of error rate if NC is unknown
N2=10; % searching limit is max(N2,nk+6)

% part 0: loading dataset
suffix = 'prefilteringAnalysis';
matfile = [suffix,'/clusterDataset.mat'];
load(matfile);

figsdir = 'clusteringResults';

% number of elements (neurons, in rows) and dimensions (measures, in columns)
[nrow, dim] = size(datanorm4Cluster);

truelabels = ones(nrow,1); % for the case it was supervised, not relevant here

% part 1: calculating euclidian distances between elements (neurons, in rows)
[dist,dmax]= similarity_euclid(datanorm4Cluster);

% part 2: pumping out the different realizations
for i=1:numberOfRepeats
    fprintf('\n>>>> Estimating the Number of Clusters, iteration=%d \n\n', i);
    classlabel = ones(nrow,maxNumberOfClusters);
    exemplars=cell(numberOfRepeats,maxNumberOfClusters);
    Rd = 'euclidean';
    % different number of clusters
    for k=minNumberOfClusters:maxNumberOfClusters
        % commented to improve speed
        %fprintf('\n>>>> running clustering algorithm with n=%d \n\n', k);
        [classlabel(:,k),centroidx] = kmeans(datanorm4Cluster,k,'replicates',numOfReplicates,'emptyaction','singleton','distance',distance);
        Q = ind2cluster(classlabel(:,k));
        ns = [];
        for j =1:numel(Q)
            ns(j) = numel(Q{j});
        end
    end
    % part 3:  estimating the number of clusters by validity indices
    SRF(i,:,:)=validity_Index_mod(datanorm4Cluster,classlabel,minNumberOfClusters,maxNumberOfClusters,truelabels,nk,N2,Rd,dist,dmax,nrow);
end

FR={'Rand','Mirkin','Hubert','Silhouette','Davies-Bouldin','Calinski-Harabasz','Hartigan','Homogeneity','Separation'};

xaxis=minNumberOfClusters:maxNumberOfClusters;

s=figure('color','none','visible','off');
hold on
set(gca,'layer','top','color','none')
for i=1:9
    ax(i)=subplot(3,3,i);
    for j=1:numberOfRepeats
      clear TT3;
      TT3=reshape(SRF(j,i,:),minNumberOfClusters,maxNumberOfClusters);
      plot(xaxis,TT3,'k','LineWidth',0.5);
      hold on;
    end
    title([FR{i}]);
    if (i>6)
      xlabel('Number of clusters');
    else
      set(gca,'XTickLabel',{})
    end
    set(gca,'fontSize',12,'LineWidth',0.5,'TickDir','out','Box','off','layer','top','color','none')
    ylimits = get(gca,'YLim');
    xrange = 5:15; % a first rough estimate of the proper range simply based on visual inspection across validation measures
    xpatch = [xrange,xrange(end:-1:1)];
    ypatch = [ylimits(1)*ones(size(xrange)),ylimits(2)*ones(size(xrange))];
    colorpatch = [.8,.8,.8];
    p=patch(xpatch,ypatch,colorpatch,'EdgeColor','none');
    uistack(p,'bottom')
end
linkaxes(ax,'x');
plot2svg([figsdir,'/numberOfClustersEstimate_',type,'.svg']);

% total elapsed time in the analysis
elapsedTime=toc(timeId);

elapsedTimeHours=floor(elapsedTime/3600);
elapsedTimeMinuts=floor(mod(elapsedTime,3600)/60);
elapsedTimeSeconds=mod(mod(elapsedTime,3600),60);

fprintf('\n>>>> total time = %d hours, %d minuts, %.3f seconds\n\n', elapsedTimeHours,elapsedTimeMinuts,elapsedTimeSeconds);
