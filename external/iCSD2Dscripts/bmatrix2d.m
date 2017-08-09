function B=bmatrix2d (nx,ny,variant)

% BMATRIX2D Calculates the matrix B.
%   Helper script for two-dimensional ICSD.

% B=bmatrix (nx,ny)
% B=bmatrix (nx,ny,variant)
%
% Calculates the matrix B for nx x ny grid.
% This matrix is used to put data from nx x ny grid
% on a larger nx+2 x ny+2 grid. 
% If optional third argument 'D' is given, 
% this is the B(D) matrix, that means the data at additional layer
% are duplicated from original data. If only two arguments are given,
% B puts zeros at the additional layer.

%   Copyright 2006-2010 Szymon Leski
%   s.leski@nencki.gov.pl

if nargin==2
    variant='B';
end;

n = nx*ny; 
M = (nx+2)*(ny+2); 

B = zeros(M,n); 

switch variant
    case 'B'
        for ii=1:M
            [xk,yk] = ind2sub([nx+2,ny+2],ii);
            if (xk==1)||(yk==1)||(xk==(nx+2))||(yk==(ny+2))
            else
                nr = sub2ind([nx ny],xk-1,yk-1);
                B(ii,nr) = 1;
            end;
        end;
    case 'D'
        for ii=1:M
            [xk,yk] = ind2sub([nx+2,ny+2],ii);
            if (xk==1)||(yk==1)||(xk==(nx+2))||(yk==(ny+2))
                if (xk==1)||(xk==(nx+2))
                    if yk==1
                        yk=2;
                    end;
                    if yk==(ny+2)
                        yk=ny+1;
                    end;
                    if (xk==1)
                        nr = sub2ind([nx ny],1,yk-1);
                    else
                        nr = sub2ind([nx ny],nx,yk-1);
                    end;
                    B(ii,nr) = 1;
                elseif (yk==1)||(yk==(ny+2))
                    if (yk==1)
                        nr = sub2ind([nx ny],xk-1,1);
                    else
                        nr = sub2ind([nx ny],xk-1,ny);
                    end;
                    B(ii,nr) = 1;
                end;
            else
                nr = sub2ind([nx ny],xk-1,yk-1);
                B(ii,nr) = 1;
            end;
        end;
    case 'E'
        for ii=1:M
            [xk,yk] = ind2sub([nx+2,ny+2],ii);
            if (xk==1)||(yk==1)||(xk==(nx+2))||(yk==(ny+2))
                if (xk==1)||(xk==(nx+2))
                    if yk==1
                        yk=2;
                    end;
                    if yk==(ny+2)
                        yk=ny+1;
                    end;
                    if (xk==1)
                        nr = sub2ind([nx ny],1,yk-1);
                    else
                        nr = sub2ind([nx ny],nx,yk-1);
                    end;
                    B(ii,nr) = 2;
                elseif (yk==1)||(yk==(ny+2))
                    if (yk==1)
                        nr = sub2ind([nx ny],xk-1,1);
                    else
                        nr = sub2ind([nx ny],xk-1,ny);
                    end;
                    B(ii,nr) = 2;
                end;
            else
                nr = sub2ind([nx ny],xk-1,yk-1);
                B(ii,nr) = 1;
            end;
        end;
end;
