function matclustsetup

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
   disp(['Compiling ',cfiles(i).name]);
   mex(cfiles(i).name)
end


disp([' ']);
disp(['You must now add ',pwd]);
disp(['to your matlab path before starting matclust.']);
disp(['Then, to start matclust, type ''matclust''']);