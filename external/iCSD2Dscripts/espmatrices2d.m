function E=espmatrices2d (nx,ny,method)

% ESPMATRICES2D Calculates E matrices for 2-D spline interpolation.
%   Helper script for two-dimensional ICSD.

% E = espmatrices2d (nx, ny, method)
% Calculates E matrices which relate value at the grid nodes
% to coefficients standing by base functions (cubic polynomials). 
% Used by initspline2d. 

%   Copyright 2006-2010 Szymon Leski
%   s.leski@nencki.gov.pl

n = nx*ny;
m = (nx-1)*(ny-1);

E = zeros(4,4,m,n);
Gx = gmatrix(nx,method);
Gy = gmatrix(ny,method);

for ii=1:m
    [xk,yk] = ind2sub([nx-1,ny-1],ii);
    nr = sub2ind([nx ny],xk,yk);
    nrx = sub2ind([nx ny],xk+1,yk);
    nry = sub2ind([nx ny],xk,yk+1);
    nrxy = sub2ind([nx ny],xk+1,yk+1);

    E(1,1, ii, nr) = 1;
    E(1,2, ii, nrx) = 1;
    E(2,1, ii, nry) = 1;
    E(2,2, ii, nrxy) = 1;
    for i=1:nx
        pk = sub2ind([nx ny],i,yk);
        E(1,3, ii, pk) = Gx(xk, i);
        E(1,4, ii, pk) = Gx(xk+1, i);
        pky = sub2ind([nx ny],i,yk+1);
        E(2,3, ii, pky) = Gx(xk, i);
        E(2,4, ii, pky) = Gx(xk+1, i);
    end;
    for i=1:ny
        pk = sub2ind([nx ny],xk,i);
        pkx = sub2ind([nx ny],xk+1,i);
        E(3,1, ii, pk) = Gy(yk, i);
        E(3,2, ii, pkx) = Gy(yk, i);
        E(4,1, ii, pk) = Gy(yk+1, i);
        E(4,2, ii, pkx) = Gy(yk+1, i);
    end;
    for i=1:ny
        for l=1:nx
            pk = sub2ind([nx ny], l, i);
            E(3,3, ii, pk) = Gx(xk, l)*Gy(yk, i);
            E(3,4, ii, pk) = Gx(xk+1, l)*Gy(yk, i);
            E(4,3, ii, pk) = Gx(xk, l)*Gy(yk+1, i);
            E(4,4, ii, pk) = Gx(xk+1, l)*Gy(yk+1, i);
        end;
    end;
end;
