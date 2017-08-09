%   Copyright 2006-2010 Szymon Leski
%   s.leski@nencki.gov.pl

function CSD = icsd2d (potentials, F, boundary)

% ICSD Calculates CSD in two dimensions. 
%
%   csd = icsd2d(potentials, F)
%   csd = icsd2d(potentials, F, boundary)
%   csd = icsd2d(potentials, Finv, 'trad')
%
%   Input: 
%
%   potentials is a 3-D array of input data
%   potentials(t,:,:) is a 2-D array of values
%   of the potential at points labelled
%   (1,1), (2,1) ... (nx,ny)
%
%   For inverse CSD:
%   F is the matrix which relates CSD to potentials.
%   (Use initstep2d, initlin2d or initspline2d to calculate F). 
%
%   boundary: 'no'='S', 'B', 'D' (optional, default = 'S')
%   Specifies the behavior of the assumed CSD at the boundary:
%   S - do nothing, B - add zeros at additional layer,
%   D - duplicate the boundary layer
%   
%   Note: If boundary is 'B' or 'D' then the matrix F
%   has to be for appropriately larger grid (nx+2 by ny+2).
%   
%   For traditional CSD:
%   Finv is the matrix calculared using inittrad2d.m, 
%   boundary must be 'trad'.
%
%
%   Output: 3-D array of CSD values at the grid points.
%   You can further interpolate it using interp2d.

%   Copyright 2006-2009 Szymon Leski
%   s.leski@nencki.gov.pl

[nt,nx,ny]=size(potentials);

if nargin==2
    boundary='S';
end;

CSD = zeros(size(potentials));
switch boundary
    case 'trad'
        CSD = FTimesData(potentials, F);
    case {'no', 'S'}
        CSD = FTimesData(potentials, inv(F));
    case 'B'
        B = bmatrix2d(nx,ny);
        R = rmatrix2d(nx,ny);
        CSD = FTimesData(potentials, inv(R*F*B));
    case 'D'
        B = bmatrix2d(nx,ny,'D');
        R = rmatrix2d(nx,ny);
        CSD = FTimesData(potentials, inv(R*F*B));
    case 'E'
        B = bmatrix2d(nx,ny,'E');
        R = rmatrix2d(nx,ny);
        CSD = FTimesData(potentials, inv(R*F*B));
end;

function out = FTimesData(fpot, invF)
[nt,nx,ny] = size(fpot);
out=permute(reshape(invF*reshape(permute(fpot,[2 3 1]),...
    [nx*ny nt]),[nx ny nt]),[3 1 2]);
% The above line is equivalent to the following code: (but faster)
% [nt,nx,ny] = size(fpot);
% n = nx*ny;
% FlatData=zeros(nt,n);     
% for k=1:n
%     FlatData(:,k)=fpot(:,k);
% end
% temp = zeros(nt,n);
% for k=1:nt
%     temp(k,:)=invF*squeeze(FlatData(k,:)');   
% end
% out = reshape(temp, size(fpot));

