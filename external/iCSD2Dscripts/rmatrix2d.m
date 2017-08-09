function R=rmatrix2d (nx,ny)

% RMATRIX2D Calculates the matrix R.
%   Helper script for two-dimensional ICSD.

% R=rmatrix2d (nx,ny)
% Calculates the matrix R for nx x ny grid.
% This matrix is used to restrict data from nx+2 x ny+2 grid
% to nx x ny grid. 

%   Copyright 2006-2010 Szymon Leski
%   s.leski@nencki.gov.pl

n = nx*ny; 
M = (nx+2)*(ny+2); 

R = zeros(n,M); 

for ii=1:n 
    [xk,yk] = ind2sub([nx,ny],ii); 
    nr = sub2ind([nx+2 ny+2],xk+1,yk+1); 
    R(ii,nr) = 1;
end;
