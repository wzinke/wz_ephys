function [rowcount] = rowCount(DATA)
%Function [rowcount] = rowCount(DATA)
%
% Function calculates how many samples (rows) DATA-matrix has in total
% (according to the first column).

% find rowcount
rowcount = length(DATA{1});