function Finv = inittrad2d(nx,ny,dx,dy, method)

% INITDELTA2D Calculates Finv matrix for 2D traditional CSD method.
%
%   Finv = inittrad2d(nx, ny, dx, dy)
%   Finv = inittrad2d(nx, ny, dx, dy, 'vaknin')
%
%   Calculates the matrix Finv which implements two-dimensional
%   Laplacian using finite-difference approximation.
%   If the third argument 'vaknin' is given, then the
%   matrix Finv implements also the Vaknin procedure for obtaining CSD at
%   the boundary. 
%
%   Input:
%
%   nx, ny - size of the grid (number of nodes)
%   dx, dy - spacing of the grid

%   Copyright 2006-2010 Szymon Leski
%   s.leski@nencki.gov.pl

if nargin==4
    method = 'trad';
end;

n = nx*ny;

Finv = zeros(n);
for ii=1:n
    [xi,yi] = ind2sub([nx,ny],ii);
    if (xi>1)&&(yi>1)&&(xi<nx)&&(yi<ny)
        iipx = sub2ind([nx ny], xi+1, yi);
        iimx = sub2ind([nx ny], xi-1, yi);
        iipy = sub2ind([nx ny], xi, yi+1);
        iimy = sub2ind([nx ny], xi, yi-1);

        Finv(ii,ii) = 2/dx^2+2/dy^2;
        Finv(ii,iipx) = -1/dx^2;
        Finv(ii,iimx) = -1/dx^2;
        Finv(ii,iipy) = -1/dy^2;
        Finv(ii,iimy) = -1/dy^2;
    end;
    if strcmpi(method, 'vaknin')
        if xi==1
            if yi==1
                iipx = sub2ind([nx ny], xi+1, yi);
                iipy = sub2ind([nx ny], xi, yi+1);
                Finv(ii,ii) = 1/dx^2 + 1/dy^2;
                Finv(ii,iipx) = -1/dx^2;
                Finv(ii,iipy) = -1/dy^2;
            elseif yi==ny
                iipx = sub2ind([nx ny], xi+1, yi);
                iimy = sub2ind([nx ny], xi, yi-1);
                Finv(ii,ii) = 1/dx^2 + 1/dy^2;
                Finv(ii,iipx) = -1/dx^2;
                Finv(ii,iimy) = -1/dy^2;
            else
                iipy = sub2ind([nx ny], xi, yi+1);
                iimy = sub2ind([nx ny], xi, yi-1);
                iipx = sub2ind([nx ny], xi+1, yi);
                Finv(ii,ii) = 1/dx^2 + 2/dy^2;
                Finv(ii,iipy) = -1/dy^2;
                Finv(ii,iimy) = -1/dy^2;
                Finv(ii,iipx) = -1/dx^2;
            end;
        end;
        if xi==nx
            if yi==1
                iimx = sub2ind([nx ny], xi-1, yi);
                iipy = sub2ind([nx ny], xi, yi+1);
                Finv(ii,ii) = 1/dx^2 + 1/dy^2;
                Finv(ii,iimx) = -1/dx^2;
                Finv(ii,iipy) = -1/dy^2;
            elseif yi==ny
                iimx = sub2ind([nx ny], xi-1, yi);
                iimy = sub2ind([nx ny], xi, yi-1);
                Finv(ii,ii) = 1/dx^2 + 1/dy^2;
                Finv(ii,iimx) = -1/dx^2;
                Finv(ii,iimy) = -1/dy^2;
            else
                iipy = sub2ind([nx ny], xi, yi+1);
                iimy = sub2ind([nx ny], xi, yi-1);
                iimx = sub2ind([nx ny], xi-1, yi);
                Finv(ii,ii) = 1/dx^2 + 2/dy^2;
                Finv(ii,iipy) = -1/dy^2;
                Finv(ii,iimy) = -1/dy^2;
                Finv(ii,iimx) = -1/dx^2;
            end;
        end;
        if (xi>1)&&(yi==1)&&(xi<nx)
            iipx = sub2ind([nx ny], xi+1, yi);
            iimx = sub2ind([nx ny], xi-1, yi);
            iipy = sub2ind([nx ny], xi, yi+1);
            Finv(ii,ii) = 2/dx^2 + 1/dy^2;
            Finv(ii,iipx) = -1/dx^2;
            Finv(ii,iimx) = -1/dx^2;
            Finv(ii,iipy) = -1/dy^2;
        end;
        if (xi>1)&&(yi==ny)&&(xi<nx)
            iipx = sub2ind([nx ny], xi+1, yi);
            iimx = sub2ind([nx ny], xi-1, yi);
            iimy = sub2ind([nx ny], xi, yi-1);
            Finv(ii,ii) = 2/dx^2 + 1/dy^2;
            Finv(ii,iipx) = -1/dx^2;
            Finv(ii,iimx) = -1/dx^2;
            Finv(ii,iimy) = -1/dy^2;
        end;
    end;
end;
