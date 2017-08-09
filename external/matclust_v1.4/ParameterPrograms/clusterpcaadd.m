function [out, names] = clusterpcaadd();

global clustdata;
global clustattrib;

dlg = inputdlg('Enter cluster number');
clustnum = str2num(dlg{1});
index = clustattrib.clusters{clustnum}.index;
cd(clustattrib.currentfilepath);
load(clustattrib.datafile);



%compute the pca's of each channel
out = [];
for i = 1:4
    covariance = cov(double(reshape(waves(:,i,index),40,[])'));
    tmpwaves = double(reshape(waves(:,i,:), 40,[]));
    coef = pcacov(covariance);
    coef = coef(:,1);
    out(:,i) = (tmpwaves'*(coef));
end


names = {'pca1','pca2','pca3','pca4'};

