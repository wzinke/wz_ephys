%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% General parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

i2mPath;
BinTab=[10  30 50 70 90 110 130]%[10 20 30];%[30 20 10 5];
SizeMax=5;
nbpatterns=300;
n_sample=1000;
N=10;

rPJ = LoadRasterPJ('PJ/PJraster1.ras');
NbNeurons=10;
f=find(rPJ(2,:)<=NbNeurons);
rPJ2 = rPJ(:,f);

for BinIndex=1:length(BinTab)
    BinSize=BinTab(BinIndex);
    
    % Load the raster data
    raster{1,1} = AlainRaster(rPJ2,BinSize,0,max(rPJ2(1,:)));
    %raster{2,1} = [];
    %raster{3,1} = [];

    % Compute the Correlation matrices and the mean from the rasters
    datalenRaw = zeros(3,1);
    for i=1:size(raster,1)
        mAv{i,BinIndex} = zeros(N,1);
        CAv{i,BinIndex} = zeros(N,N);
        C1Av{i,BinIndex} = zeros(N,N);
        for j=1:size(raster,2)
            if ~isempty(raster{i,j})
                datalenRaw(i,j) = size(raster{i,j},2);
                mRaw{i,j} = mean(raster{i,j},2);
                CRaw{i,j} = raster{i,j}*transpose(raster{i,j})/datalenRaw(i,j);
                C1Raw{i,j} = raster{i,j}(:,1:(datalenRaw(i,j)-1))*transpose(raster{i,j}(:,2:datalenRaw(i,j)))/(datalenRaw(i,j)-1);                            
                mAv{i,BinIndex} = mAv{i,BinIndex} + datalenRaw(i,j)*mRaw{i,j};
                CAv{i,BinIndex} = CAv{i,BinIndex} + datalenRaw(i,j)*CRaw{i,j};
                C1Av{i,BinIndex} = C1Av{i,BinIndex} + datalenRaw(i,j)*C1Raw{i,j};
            end
        end
        datalen(i) = sum(datalenRaw(i,:));
        mAv{i,BinIndex} = mAv{i,BinIndex}/datalen(i);
        CAv{i,BinIndex} = CAv{i,BinIndex}/datalen(i);
        C1Av{i,BinIndex} = C1Av{i,BinIndex}/datalen(i);
    end
    MatCorr{BinIndex} = C1Av{i,BinIndex}./CAv{i,BinIndex};
    CorrIndex(BinIndex) = mean(mean(C1Av{i,BinIndex}./CAv{i,BinIndex}));
end