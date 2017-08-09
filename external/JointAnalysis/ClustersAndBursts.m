
BurstSp=[];
SeparateSp=[];

for l=1:length(CellsInterested)
    ChannelToAnalyze = (CellsInterested(l)); %write here the order of the channel in DataCell
    
    z=unique(cluster_class(1:end,1));
    ClustersInBursts=zeros(1,length(z));
    ClustersOutBursts=zeros(1,length(z));
    
    a = min(find(channels==ChannelToAnalyze));%"channels" is calculated earlier by "DataCellSpikeSorter.m"
    b = max(find(channels==ChannelToAnalyze));
    
    SpikeTypes=TypeOfSpikeCell{l};    
    
    for j = a : b
        if SpikeTypes(j-a+1)==1
            if sum(ismember(z,0))~=0 %check if there is unclustered "0" type spike, if there is start index from 0+1
                ClustersInBursts(1,cluster_class(j,1)+1) = ClustersInBursts(1,cluster_class(j,1)+1)+1;
            else
                ClustersInBursts(1,cluster_class(j,1)) = ClustersInBursts(1,cluster_class(j,1))+1;
            end
        else
            if sum(ismember(z,0))~=0 %check if there is unclustered "0" type spike, if there is start index from 0+1
                ClustersOutBursts(1,cluster_class(j,1)+1) = ClustersOutBursts(1,cluster_class(j,1)+1)+1;
            else
                ClustersOutBursts(1,cluster_class(j,1)) = ClustersOutBursts(1,cluster_class(j,1))+1;
            end
        end
        
    end
    
    BurstSp =[BurstSp;ClustersInBursts];
    SeparateSp=[SeparateSp;ClustersOutBursts];
end
