function output = clustExtract(filename)

% output = clustExtract(filename)
%
% this funtion returns a structure array with fields timestamps and
% parameterindex for each cluster in the matclust file.

load(filename);
output = [];
for i = 1:length(clustattrib.clusters)
    output(i).timestamps = clustdata.params(clustattrib.clusters{i}.index,1);
    output(i).parameterindex = clustattrib.clusters{i}.index;
end