function [DATA] = medianFilterData(DATA, winlen, column)
%Function [DATA] = medianFilterData(DATA, winlen, column)
%
% Performs median filtering to the datapoints in the specified column on
% DATA-structure. Further information see help medianFilter.
% Column specifies the column to filter the data with

DATA{column} = medianFilter(DATA{column}, winlen);