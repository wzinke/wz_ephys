function [out, names] = testparamadd();

global clustdata;

out = clustdata.params(:,1:2);
names = {'this','that'};
