function [out, names] = pcaadd();

global clustdata;
global clustattrib;

dlg = inputdlg('Enter threshold (uV)');
thresh = str2num(dlg{1});
cd(clustattrib.currentfilepath);
load(clustattrib.datafile);

amps = clustdata.params(:,2:5);
pass = find( (amps(:,1) >= thresh) | (amps(:,2) >= thresh) | (amps(:,3) >= thresh) | (amps(:,4) >= thresh) );
out = spikeprincecomp(waves, 6, amps, pass);
names = {'pca1','pca2','pca3','pca4','pca5','pca6'};

