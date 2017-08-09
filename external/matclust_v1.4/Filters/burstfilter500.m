function out = NonComplexSpikes()

global clustdata;
global clustattrib;

clustnum = inputdlg('Which cluster?','Filter non-complex spikes');
clustnum = str2num(clustnum{1});
historywindow = 1000; 
if ~isempty(clustnum)
  
    out = false(length(clustdata.params(:,1)),1);
    index = clustattrib.clusters{clustnum}.index;
    [ma, maxchannel] = max(var(clustdata.params(index,2:5)));
    
    if length(index) > 2
        timediff = diff(clustdata.params(index,1));
        inpoints = find((timediff> 5000)&(timediff>5000))+1;
        out(index(inpoints)) = true;
    end
    
    
else
    out = [];
end
    
    
    