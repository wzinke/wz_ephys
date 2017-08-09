function labels=cluster2label(clusters,numElements)

numClusters=length(clusters);
labels=zeros(numElements,1);
for i=1:numClusters
  labels(clusters{i}(:))=i;
end