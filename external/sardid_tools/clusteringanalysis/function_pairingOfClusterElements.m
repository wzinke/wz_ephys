function labels = function_pairingOfClusterElements(type,datanorm4Cluster)

% This function runs different realizations of the Kmeans/Kmedians clustering
% algorithm with a previously established number of clusters to evaluate
% to what extend clusters are composed of the same elements

% type = 'Kmeans';
% type = 'Kmedians';
if strcmp(type,'Kmeans')
  distance = 'sqEuclidean';
elseif strcmp(type,'Kmedians')
  distance = 'cityblock';
else
  error('type of clustering has to be either Kmeans or Kmedians');
end

% parameters
numberOfRepeats=500;      % number of realizations
numOfReplicates=50;       % number of replicates

% part 0: loading dataset (to modify the dataset, please go to measuresPrefiltering.m)

% REMEMBER to modify the suffix accordingly!
% suffix = 'WF_ACT_wholeTrial';
% matfile = ['../',suffix,'/clusterDataset_',suffix,'.mat'];
% load(matfile);

% number of clusters to consider (based on estimateNumberOfClusters)
% numC=5:10; % 5:15 % no need to go that further, as cells are not consistently grouped together for >7 clusters
numC = 7; % this is only to validate results for best num of clusters (7)
cutoff_clusterElements = 5; % at least 5 elements in a cluster to be considered

% number of elements (neurons, in rows) and dimensions (measures, in columns)
[nrow, dim] = size(datanorm4Cluster);

% initializing labels (clusters to which neurons belong to)
labels=zeros(length(numC),nrow);

ind=0;
for numClusters=numC
    ind=ind+1;
    fprintf('\n>>>> starting clustering algorithm with %d clusters \n\n', numClusters);

    pairwiseCluster=tril(nan(nrow,nrow),0)+triu(zeros(nrow),1); % initializing the matrix of pairwise same cluster belonging

    classlabel = zeros(nrow,numberOfRepeats); % it keeps record of neurons' cluster through realizations

    numOfIntersections=zeros(numClusters,numClusters);
    minNumOfIntersections=nrow*ones(numClusters,1);
    maxNumOfElem=zeros(numClusters,1);

    for i=1:numberOfRepeats

        % part 1: running clustering algorithm

        if (mod(i,50)==0)
            fprintf('\n++++ running clustering algorithm, iteration=%d \n\n', i);
        end

        [classlabel(:,i),centroidx] = kmeans(datanorm4Cluster,numClusters,'replicates',numOfReplicates,'emptyaction','singleton','distance',distance);

        if(i>1)
          Qprev=Q;
        end
        Q = ind2cluster(classlabel(:,i));

        % part 2: control for neurons belonging to same clusters

        ns = [];
        for j =1:numel(Q)
            ns(j) = numel(Q{j});
            for k=1:ns(j)
                for l=k+1:ns(j)
                    pairwiseCluster(Q{j}(k),Q{j}(l))=nansum(pairwiseCluster(Q{j}(k),Q{j}(l))+1);
                end
            end
        end
    end

    pairwiseCluster=pairwiseCluster/numberOfRepeats;
    symPairwiseCluster=tril(pairwiseCluster')+triu(pairwiseCluster);

    cutoff=0.9;
    symPairwiseClusterMod=zeros(size(pairwiseCluster));
    symPairwiseClusterMod(symPairwiseCluster>=cutoff)=symPairwiseCluster(symPairwiseCluster>=cutoff);
    reorderedSymPairwiseCluster=symrcm(symPairwiseClusterMod);

%      if ind == 1
%        s=figure('color','none','visible','off');
%        hold on
%        set(gca,'layer','top','color','none')
%        [cc1,hc1] =contourf(symPairwiseClusterMod);
%        set(hc1,'LineStyle','none');
%        colormap copper
%        set(gca,'fontSize',16,'LineWidth',1,'TickDir','out','Box','off')
%        plot2svg([suffix,'/pairwise',num2str(numClusters),'clusterSparse_',type,'_final.svg']);
%      end
%
%      s=figure('color','none','visible','off');
%      hold on
%      set(gca,'layer','top','color','none')
%      [cc2,hc2] = contourf(symPairwiseClusterMod(reorderedSymPairwiseCluster,reorderedSymPairwiseCluster));
%      set(hc2,'LineStyle','none');
%      colormap copper
%      set(gca,'fontSize',16,'LineWidth',1,'TickDir','out','Box','off')
%      plot2svg([suffix,'/pairwise',num2str(numClusters),'clusterReordered_',type,'_final.svg']);

    pairwiseClusterReordered=triu(symPairwiseClusterMod(reorderedSymPairwiseCluster,reorderedSymPairwiseCluster));

    k=0;
    allElements=[];
    for i=1:nrow
        rowElements=find(pairwiseClusterReordered(i,:)>=cutoff);
        if(~isempty(rowElements))
            rowElements=union(i,rowElements);
            if(isempty(intersect(allElements,rowElements)))
                k=k+1;
                clust{k}=rowElements;
            else
                clust{k}=union(clust{k},rowElements);
            end
            allElements=union(allElements,rowElements);
        end
    end

    % check again that all clusters don't share elements (sometimes they do, see error below)
    flag=1;
    while (flag==1)
        flag=0;
        l=length(clust);
        for i=1:(length(clust)-1)
            if(isempty(clust{i}))
                clust{i}=clust{i+1};
                clust{i+1}={};
            elseif(~isempty(intersect(clust{i},clust{i+1})))
                clust{i}=union(clust{i},clust{i+1});
                clust{i+1}={};
                l=l-1;
                flag=1;
            end
        end
        clust=clust(1:l);
    end

    cont=0;
    allElements=[];
    for i=1:length(clust)
        if(numel(clust{i})>=cutoff_clusterElements)
            cont=cont+1;
            for j=1:numel(clust{i})
                clustFilt{ind,cont}(j)=reorderedSymPairwiseCluster(clust{i}(j));
                allElements=[allElements,clustFilt{ind,cont}(j)];
            end
            clustFilt{ind,cont}=sort(clustFilt{ind,cont});
        end
    end
    clear clust;

    percElem=100*length(allElements)/nrow;
    display(percElem);
    percElements(ind)=100*length(allElements)/nrow;

    labels(ind,:)=cluster2label(clustFilt(ind,:),nrow);
end
