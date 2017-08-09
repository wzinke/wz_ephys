%% Preparing format of the data for further spike sorting with Wave_Clus
% 
%  Data is stored in previously formed DataCell (a self created data format
%  for our lab to store MEA data). DataCell has seven (with the raw data)or
%  six columns (only with spike times) according to user's preference. 
%  We use in this code the second column where spike times of each channel
%  recorded in milliseconds and the sixth column where the spike waveforms 
%  belong to those channels are stored.

numSpikes = sum(cell2mat(DataCell(1:1,2:2))); % from 9x6 = 54 channels for 6-well MEA
WF_length = size (DataCell{1,6},2); %default WF_length is 64
index = zeros(1,numSpikes); % spike times
spikes = zeros(numSpikes,WF_length); % spike waveforms
channels = zeros(1,numSpikes); % channels where corresponding spike WFs occurs

k=0;

for i = 1 : size(DataCell,1)
    
    numberOFspikes = length(DataCell{i,3});
    index (k+1:k+numberOFspikes) = DataCell{i,3};
    spikes (k+1:k+numberOFspikes, 1:WF_length) = DataCell{i,6};
    channels (k+1:k+numberOFspikes) = i;
    k=k+numberOFspikes;
    
end


ChannelNames = cell(1,size(DataCell,1));
for k = 1 : size(DataCell,1)
ChannelNames{k} = DataCell{k,1}(end-1:end);
end

