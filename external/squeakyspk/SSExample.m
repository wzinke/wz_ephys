% Example script for usage of SqueakySpk in data processing and storage.
% I suggest running this script cell by cell to get an idea of the workflow
% 
% Jon Newman
% 2010-06-23

% Object FID - Name for the SS object
fid = '20100604_16269_SS';

% Load data into data structs
spkdat = load('testdata.mat');

% Examine waveforms before processing. Look at all the artifacts.
figure()
plot(spkdat.waveform);
xlabel('Time (samp)')
ylabel('mV')

%% Instantiate a SqueakySpk Object
SStest = SqueakySpk(fid,spkdat,25000);

% Remove meaningless channels
SStest.RemoveChannel();

% Remove spikes with blanks
SStest.RemoveSpkWithBlank();

% Perform a hard p2p cut at 175 uV
SStest.HardThreshold(175);

% Clustering
SStest.WaveClus(3,20,'wav',1);
SStest.RemoveUnit(0); % remove unsorted data

% Supervised unit deletion by average waveform
SStest.WeedUnitByWaveform()

%% Examime some data to make sure results of sorting and cleaning look good
SStest.RasterWave_Comp([],'both');

% Return clean data to workspace
cleandata = SStest.ReturnClean();

% Examine the improvement after processing
figure()
plot(cleandata.cwaveform);
xlabel('Time (samp)')
ylabel('mV')
%% Save the SStest data object to your working directory
save(fid, 'SStest');