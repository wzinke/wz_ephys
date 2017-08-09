function ExcludeOverlappingPoints(clustnum)

global clustdata;
global clustattrib;
global figattrib;

overlap = analyzeoverlap(0);
clustindex = find((overlap.list(:,1) == clustnum)|(overlap.list(:,2) == clustnum));
for i = 1:length(clustindex)
    currindex = overlap.points{i};
    tmpadd = currindex(:);
    tmpadd(:,2) = 0;
    tmpadd = int32(tmpadd);
    tmpadd(:,2) = fastbitset(tmpadd(:,2),2,true);

    try
        clustattrib.eventeditindex{clustnum} = [clustattrib.eventeditindex{clustnum};tmpadd];
    catch
        clustattrib.eventeditindex{clustnum} = tmpadd;
    end
    
end

in = fastorbit(clustattrib.pointexclude,clustattrib.clusters{clustnum}.polyindex(:,4)); %which points are excluded at any polygon
in2 = fastorbit(clustattrib.pointinclude,clustattrib.clusters{clustnum}.polyindex(:,4)); %which points have been included at any polygon
in3 = false(length(clustdata.params(:,1)),1); %which points are excluded individually

try
    if ~isempty(clustattrib.eventeditindex{clustnum})
        in3(clustattrib.eventeditindex{clustnum}(:,1)) = true; %if any points were excluded individually, it will be stored in clustattrib.eventeditindex
    end
end

clustattrib.clusters{clustnum}.index = uint32(find(~in & in2 & ~in3));

matclust('addnewstate','exclude overlap',figattrib.handles); 
matclust('plotgraph',figattrib.handles); 


