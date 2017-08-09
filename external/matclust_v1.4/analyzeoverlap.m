function overlap = analyzeoverlap(displaymode)

global clustattrib;

if (displaymode)
   disp(' ')
   disp('Overlaps:')
   disp(' ')
end
overlap = [];
count = 1;
for i = 1:length(clustattrib.clusters)
      index = [];
      try   
         index = clustattrib.clusters{i}.index;   
      end
      for j = i+1:length(clustattrib.clusters)
         if (j~=i)
            index2 = [];
            try
               index2 = clustattrib.clusters{j}.index;   
            end
            opoints = intersect(index, index2);
            if ~isempty(opoints)
               overlap.list(count,1:3) = [i j length(opoints)];
               overlap.points{count} = opoints;
               if (displaymode)
                  disp(['Cluster ',num2str(i),' overlaps with cluster ',num2str(j),': ',num2str(length(opoints)),' points.'])
               end
               count = count+1;
            end
         end
      end
end