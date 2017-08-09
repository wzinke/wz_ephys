
clear;
close all;
clc;

type = 'Kmeans';
resultsdir = 'clusteringResults';
matfile = [resultsdir,'/',type,'ClusteringResults'];
load(matfile);

[numNeurons,numMeasures] = size(datanorm4Cluster);

figure('color','none','visible','off')
hold on
set(gca,'layer','top','color','none')
p=mypcolor(1:numMeasures,1:numNeurons,datanorm4Cluster);
rgb_red=[1,0.35,0.35]; % red
rgb_blue=[0.35,0.35,1]; % blue
map = diverging_map(0:1/255:1,rgb_red,rgb_blue);
colormap(map)
set(p,'EdgeColor','interp');
axis tight
axis off
filename=[resultsdir,'/rawHeatMap'];
svgFile=[filename,'.svg'];
plot2svg(svgFile);

% resorting data according to the dendrogram of the k-means clustering
numClusters=7;
offset=4; % it starts with numClusters=5
k=1;
boundaries=[];
sortedData=zeros(size(datanorm4Cluster));
numCellType = nan(size(celltype));
numCellType(strcmp('narrow',celltype)) = 1;
numCellType(strcmp('broad',celltype)) = 2;
numCellType(strcmp('fuzzy',celltype)) = 1.5;
sortedCelltype=zeros(size(numCellType));
for i=clustOrder
    if(~isempty([clustFilt{numClusters-offset,i}]))
        tmp = datanorm4Cluster([clustFilt{numClusters-offset,i}],:);
        [tmp2,ind]=sort(numCellType([clustFilt{numClusters-offset,i}]));
        sortedCellTypes(k:k+length([clustFilt{numClusters-offset,i}])-1)=tmp2;
        sortedData(k:k+length([clustFilt{numClusters-offset,i}])-1,:)=tmp(ind,:);
        k=k+length([clustFilt{numClusters-offset,i}]);
        boundaries=[boundaries;k];
    end
end
boundaries(end)=[];

figure('color','none','visible','off')
hold on
set(gca,'layer','top','color','none')
p=mypcolor(1:numMeasures,1:numNeurons,sortedData);
colormap(map);
colorbar
set(p,'EdgeColor','interp');
linesx=ones(length(clustOrder)-1,1)*[1 numMeasures+1];
linesy=[boundaries boundaries];
plot(linesx',linesy','k:','LineWidth',1)
axis tight
axis off
filename=[resultsdir,'/colorMap'];
svgFile=[filename,'.svg'];
plot2svg(svgFile);

figure('color','none','visible','off')
hold on
set(gca,'layer','top','color','none')
p=mypcolor(1:numMeasures,1:numNeurons,sortedData);
colormap(map);
set(p,'EdgeColor','interp');
linesx=ones(length(clustOrder)-1,1)*[1 numMeasures+1];
linesy=[boundaries boundaries];
plot(linesx',linesy','k:','LineWidth',1)
axis tight
axis off
filename=[resultsdir,'/sortedHeatMap'];
svgFile=[filename,'.svg'];
plot2svg(svgFile);

figure('color','none','visible','off')
hexcolors={'1D0091','007BFF','03FC5E','FFFF00','FFAE00','FF4000','000000','A8100D'};
colororder=zeros(length(hexcolors),3);
for i=1:length(hexcolors)
    colororder(i,:) = rgbconv(hexcolors{i});
end
set(gca, 'ColorOrder', colororder);
hold on
set(gca,'layer','top','color','none')
linesx=ones(length(clustOrder),2);
linesy=[1 boundaries' ;boundaries' length(numCellType)]';
plot(linesx',linesy','LineWidth',5)
axis tight
axis off
filename=[resultsdir,'/clusterColors'];
svgFile=[filename,'.svg'];
plot2svg(svgFile);

figure('color','none','visible','off')
hold on
set(gca,'layer','top','color','none')
p=mypcolor(1,1:length(numCellType),sortedCellTypes');
map=[
  1 0 0
  .5 .5 .5
  0 0 1
];
colormap(map);
axis tight
xlim([0 7.5]);
axis off
filename=[resultsdir,'/cellTypeColors'];
svgFile=[filename,'.svg'];
plot2svg(svgFile);

save(matfile,'numNeurons','numMeasures','numCellType','clustOrder','hexcolors','-append');
