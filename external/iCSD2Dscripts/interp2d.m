function out = interp2d (csd, VX, VY, method, boundary)

% INTERP2D Interpolates two-dimensional CSD data.
%
%   out = interp2d(csd, VX, VY, method)
%   out = interp2d(csd, VX, VY, method, boundary)
%
%   Interpolates two-dimensional CSD data between nodes 
%   taking care of the behavior at the boundary.
%
%   Input:
%
%   csd - 3-D array of CSD at the nodes, 
%   for each t csd(t,:,:) is a 2-D matrix of CSD values 
%   at points labelled
%   (1,1), (2,1), ... (nx, ny)
%
%   VX, VY - these vectors define points at which
%   the interpolated values are calculated
%   (the values in VX range from 1 to nx, in VY - from 1 to ny).
%
%   method - 'step', 'lin', 'splinen' or 'splinem'
%   sets the interpolation method (nearest neighbor, linear,
%   natural spline, not-a-knot spline)
%
%   boundary: 'no'='S', 'B', 'D' (optional, default = 'S')
%   Specifies the behavior of the assumed CSD at the boundary:
%   S - do nothing, B - add zeros at additional layer,
%   D - duplicate the boundary layer
%
%   Output:
%   3-D array of size: size(csd,1) x length(VX) x length(VY)
%   containing interpolated values of CSD.
%
%   NOTE: natural spline interpolation is slow in present version.

%   Copyright 2006-2010 Szymon Leski
%   s.leski@nencki.gov.pl

if nargin==4
    boundary='S';
end;

[nt,nx,ny] = size(csd);
out = zeros(nt, length(VX), length(VY));

if strcmp(method(1:2), 'st')
    for l = 1:nt
        out(l,:,:,:) = interp2(squeeze(csd(l,:,:)), VY, VX', 'nearest');
    end;
elseif strcmp(method(1:3), 'lin')
    for l = 1:nt
        out(l,:,:,:) = interp2(squeeze(csd(l,:,:)), VY, VX', 'linear');
    end;
elseif strcmp(method, 'splinen')
    switch boundary
        case {'no', 'S'}
            E=espmatrices2d(nx,ny,'n');
        case {'B', 'D', 'E'}
            E=espmatrices2d(nx+2,ny+2,'n');
    end;
    for l = 1:nt
        out(l,:,:) = nsplint2d(squeeze(csd(l,:,:)), VX, VY, E, boundary);
    end;
elseif strcmp(method, 'splinem')
    switch boundary
        case {'no', 'S'}
            for l = 1:nt
                out(l,:,:,:) = interp2(squeeze(csd(l,:,:,:)),VY, VX', 'spline');
            end;
        case {'B', 'D', 'E'}
            if boundary=='B'
                B = bmatrix2d(nx, ny);
            end;
            if boundary=='D'
                B = bmatrix2d(nx, ny, 'D');
            end;
            if boundary=='E'
                B = bmatrix2d(nx, ny, 'E');
            end;
            for l = 1:nt
                out(l,:,:,:) = interp2(0:ny+1, 0:nx+1, ...
                reshape(B*csd(l,:)',nx+2,ny+2), VY, VX', 'spline');
            end;
    end;
end;
