function [DATA] = concatenateData(data1, data2)
% function [DATA] = appendData(data1, data2)
%  Appends datapoints from data2 after the last datapoint of data1. 
%  Use carefully, because the continuity of your data is in danger! 

disp('Concatenating two data chunks...')

rowcount = rowCount(data1);
disp(['Datamatrix contains ' num2str(rowcount) ' rows before operation (1.st column).']);

for i=1:size(data1, 2)
   DATA{i} = [data1{i};data2{i}];
end

rowcount = rowCount(DATA);
disp(['Datamatrix contains ' num2str(rowcount) ' rows after operation (1.st column).']);