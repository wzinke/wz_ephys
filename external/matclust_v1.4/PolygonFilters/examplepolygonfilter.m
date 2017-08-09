function result = examplepolygonfilter(vertices,ax,data)
%result = examplepolygonfilter(vertices,axes,data)
%An example filter for matclust polygon filter functions
%VERTICES, [X Y], is an n by 2 matrix, where n is the number of vertices in the input polygon
%AX, [a1 a2], are the axes indices into DATA representing the x and y dimensions in which the polygon was drawn 
%DATA, is an n by m matrix where each row is a point and each column is a dimension


insidePolygon = fastinpoly(data(:,ax(1)),data(:,ax(2)),vertices(:,1),vertices(:,2)); %find the points of the current filter that are inside the polygon

meanpoint = mean(data(find(insidePolygon),2:end)); %find the center of the polygon in the full space, excluding the fist column which is usually time
distRank = pdist2(meanpoint, data(:,2:end))'; %for each point, calculate the euclidean distance to the center
result = (distRank < prctile(distRank, 25)); %keep points within the 25th percentile for distance


