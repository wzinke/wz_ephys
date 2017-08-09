function [out, names] = pcaadd();

global clustdata;
global clustattrib;

dlg = inputdlg('Enter threshold (uV)');
thresh = str2num(dlg{1});
cd(clustattrib.currentfilepath);
load(clustattrib.datafile);

amps = clustdata.params(:,2:5);
pass = find( (amps(:,1) >= thresh) | (amps(:,2) >= thresh) | (amps(:,3) >= thresh) | (amps(:,4) >= thresh) );

%compute the pca's of each channel
out = [];
for i = 1:4
    covariance = cov(double(reshape(waves(:,i,pass),40,[])'));
    tmpwaves = double(reshape(waves(:,i,:), 40,[]));
    coef = pcacov(covariance);
    coef = coef(:,2);
    out(:,i) = (tmpwaves'*(coef));
end

names = {'pca1','pca2','pca3','pca4'};

