function nonnetworksetup

currdir = pwd;

tmploc = strfind(which('setuplocal'),'.m')-11;
filelocation = which('setuplocal');
filelocation = filelocation(1:tmploc);  %contains the path to MatClust
cd(filelocation);



load defaultsbin;

matclust_defaults.Cluster0Color = defaults.Cluster0Color;
matclust_defaults.ClusterColors = defaults.ClusterColors;
matclust_defaults.GraphBackgroundColor = defaults.GraphBackgroundColor;
matclust_defaults.MaxUndos = defaults.MaxUndos;
matclust_defaults.ResFactor = defaults.ResFactor;
matclust_defaults.UnitsPerSec = defaults.UnitsPerSec;
matclust_defaults.DataFolder = pwd;
matclust_defaults.UserFolder = pwd;

try
 save matclust_defaults matclust_defaults -v6; 
catch
 save matclust_defaults matclust_defaults;
end

cfiles = dir('*.c');
for i = 1:length(cfiles)
   mex(cfiles(i).name)
end