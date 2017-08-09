function out = noisefilter()
% This filter allows only the points where at least one one the amplitude
% channels went above the set threshhold to pass through the filter.  If the threshhold was set too
% low during the experiment, this filter can speed up the plotting function
% for matclust drastically.  
global clustdata;

a = inputdlg('Enter thresh (uV)');
thresh = str2num(a{1});

out = ((clustdata.params(:,2)>thresh)|(clustdata.params(:,3)>thresh)|(clustdata.params(:,4)>thresh)|(clustdata.params(:,5)>thresh));




    
    
    