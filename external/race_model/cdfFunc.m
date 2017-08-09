function [cdfarr, edges] = cdfFunc(arr, binwidth, min, max);
%% cdfFunc.m
% cumulative distribution function

edges = [min:binwidth:max];

HistArr = histc(arr, edges);

x=1; st = 1; cdfarr = []; proparr = [];
for x = 1:1:length(HistArr)
    cdfarr(st) = sum(HistArr(1:x))/sum(HistArr);
    proparr(st) = HistArr(x)/sum(HistArr);
    st = st + 1;
end