function c = clusterStartScript()

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;
handles = figattrib.handles;

currdir = pwd;
cd ..
path = pwd;
cd(currdir);
matclust('importtimefilt',path,'times.mat');
matclust('setZoomRange',99.9);
matclust('zoomAllOut');
