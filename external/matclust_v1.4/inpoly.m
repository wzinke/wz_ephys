function in = inpoly(x,y,polyX,polyY)

%if ((polyY(1) == polyX(end))&(polyY(1) == polyY(end)))
    polyX = polyX(1:end-1);
    polyY = polyY(1:end-1);
%end
triangles = [];
tri = triangulate(polyX,polyY);
tri = tri';
n = size(tri,2);
tri = reshape(tri,[3 1 n]);
tri = [tri tri];
tri(:,1,:) = polyX(tri(:,1,:));
tri(:,2,:) = polyY(tri(:,2,:));
in = intriangle([x y],tri);