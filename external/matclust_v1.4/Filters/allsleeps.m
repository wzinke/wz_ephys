function out = allsleeps()

global clustdata;
global clustattrib;


names = lower(clustdata.timefilternames);
hits = strfind(names, 'sleep');
hitvector = [];
out = false(length(clustdata.params(:,1)),1);
for i = 1:length(hits);
   if (~isempty(hits{i}))
      hitvector(i) = 1;
   else
      hitvector(i) = 0;
   end
end


hitindex = find(hitvector);
for i = hitindex
   out(find((clustdata.params(:,1) >= clustdata.timefilterranges(i,1)) & (clustdata.params(:,1) <= clustdata.timefilterranges(i,2)))) = true;
end



    
    
    