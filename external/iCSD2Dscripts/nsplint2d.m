function out = nsplint2d(dane, VX, VY, E, boundary)

% NSPLINT2D Natural spline interpolation in 2-D. 
%   Helper script for two-dimensional ICSD.

% VI = nsplint(data, VX, VY, E, boundary)
% Used by interp2d.
% Input: data - 2-D matrix of function values at points
% (1,1), (2,1), ... (nx, ny)
% VX, VY - these vectors define points at which
% the interpolated values are calculated
% E - the matrix E 
% boundary:
% 'no' or 'S'  - plain interploation using 'data' matrix
% B - interpolation using 'data' matrix with additional layer of zeros
% D - as B but with duplicated values at the additional layer
% For B and D the E matrix should be for the larger grid
% (eg. 6x7 instead of 4x5)
%
% NOTE: this is very simple and very slow implementation. 

%   Copyright 2006-2010 Szymon Leski
%   s.leski@nencki.gov.pl

out = zeros(length(VX),length(VY));
[nx,ny] = size(dane);
A=zeros(4,4,(nx-1)*(ny-1));
if strcmp(boundary, 'no')||strcmp(boundary, 'S')
    for i=1:4
        for j=1:4
            A(i,j,:)= squeeze(E(i,j,:,:))*dane(:);
        end;
    end;
else
    if strcmp(boundary, 'B')
        B = bmatrix2d(nx,ny);
    end;
    if strcmp(boundary, 'D')
        B = bmatrix2d(nx,ny,'D');
    end;
    RA = rmatrix2d(nx-1, ny-1);
    for i=1:4
        for j=1:4
            A(i,j,:)= RA*squeeze(E(i,j,:,:))*B*dane(:);
        end;
    end;
end;

for a = 1:length(VX)
    for b = 1:length(VY)
        out(a,b) = inpol(A, VX(a), VY(b), nx,ny);
    end;
end;

function value = inpol(A, x,y, nx,ny)

bx=0;
by=0;
if x>=nx 
    x = x-1;
    bx = 1;
end;
if y>=ny 
    y = y-1;
    by = 1;
end;
nr = sub2ind([nx-1 ny-1],floor(x),floor(y)); 

dx = x-floor(x)+bx; 
dy = y-floor(y)+by;

XX(1) = 1-dx;           
XX(2) = dx;             
XX(3) = (XX(1)^3-XX(1))/6;    
XX(4)= (XX(2)^3-XX(2))/6;     
YY(1) = 1-dy;
YY(2) = dy;
YY(3) = (YY(1)^3-YY(1))/6;
YY(4) = (YY(2)^3-YY(2))/6;

value = 0;

for i=1:4
    for j=1:4
        value = value + A(j,i,nr)*XX(i)*YY(j);
    end;
end;
