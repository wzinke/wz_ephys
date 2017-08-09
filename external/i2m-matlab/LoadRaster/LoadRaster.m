function rast = LoadRaster(namefile)
%Load a text file where the spike times have the following format:
%[time 1 of cell 1]
%[time 2 cell 1]
% ...
%[time n cell 1]
%[time 1 cell 2]
%[time 2 cell 2]
%...
% We assume here that the times are in seconds, and convert them to ms
% This is, most of the time, easy to convert the data to that format

%output:
%rast(1,k): the time of the k-th spike
%rast(2,k): the cell index of the k-th spike
times=load(namefile,'ascii');

rast(1,:) = times*1000;%Conversion to milliseconds

j=1;
rast(2,1) = j;

for i=2:length(times)
    if times(i) < times(i-1)
        j = j+1;
    end
    rast(2,i) = j;
end
