function initlin2d(name,nx,ny,dx,dy,h,zprofile,dsp)

% INITLIN2D Calculates F matrix for 2-D linear ICSD method.
%
%   initlin2d(name, nx, ny, dx, dy, h)
%   initlin2d(name, nx, ny, dx, dy, h, dsp)
%   initlin2d(name, nx, ny, dx, dy, h, zprofile)
%   initlin2d(name, nx, ny, dx, dy, h, zprofile, dsp)
%
%   Calculates the matrix F which relates CSD to potentials.
%   F(k,l) is the potential generated at k-th node by 
%   unit CSD placed at l-th node of 2-D grid if linear
%   interpolation of CSD between nodes is assumed.
%   The nodes (x,y) are numbered as follows: 1) (1,1) 
%   2) (2,1) ... n) (nx,ny). 
%
%   Input:
%
%   name is the name under which the result will be saved
%   (in .mat file in subdirectory data)
%
%   nx, ny - size of the grid (number of nodes)
%   dx, dy - spacing of the grid
%
%   The argument h controls the assumed thickness of CSD. 
%   By default it is (half) the size of the step function 
%   in z direction. If optional argument zprofile='gauss' is given
%   then the step function is replaced with a Gaussian and h is the
%   standard deviation.
%
%   The argument dsp is optional and describes the
%   geometry of electrode (potential) array if it is not rectangular
%   dsp(1..nx, 1..ny, 2)
%   dsp(:,:,1) is x-displacement
%   dsp(:,:,2) is y-displacement
%   Note that CSD is still spanned on the rectangular grid.
%
%   Output:
%   a .mat file in subdirectory data, containing F and some other variables

%   Copyright 2006-2010 Szymon Leski
%   s.leski@nencki.gov.pl

datafolder = 'data';

n = nx*ny;

parseinpargs

E = EMatrix(nx,ny);

disp(' Generating F matrix for linear distribution...');
if rectgrid
    Ft = FdMatrix_rect(nx, ny, dx, dy, h, zprofile); 
    % faster, rectangular only
else
    Ft = FdMatrix(nx, ny, dx, dy, h, zprofile, dsp); 
    % slower, general geometry
end;

F = zeros(n);
for ii=1:4
    F = F + squeeze(Ft(ii,:,:))*squeeze(E(ii,:,:));
end;

if rectgrid
    save(fname, 'F', 'name', 'dx', 'dy', 'h', ...
        'nx', 'ny', 'zprofile');
else
    save(fname, 'F', 'name', 'dx', 'dy', 'h', ...
        'nx', 'ny', 'dsp', 'zprofile');
end;

function Fd=FdMatrix (nx,ny, dx, dy, h, zprofile, dsp)
n = nx*ny;
m = (nx-1)*(ny-1);
Fd = zeros(4,n,m);
counter_1 = 0;
for jj=1:m
    for ii=1:n
        [xi,yi] = ind2sub([nx,ny],ii);
        [xj,yj] = ind2sub([nx-1,ny-1],jj);
        xt = (xi - xj)*dx-dsp(xi,yi,1); % For non-rectangular 
        yt = (yi - yj)*dy-dsp(xi,yi,2); % grid of electrodes
        switch zprofile
            case 'step'
                for a = 1:4
                    Fd(a,ii,jj) = 1/2/pi*dblquad(@(x,y) ...
                        fun(x,y,xt,yt, a, dx, dy, h), 0,dx,0,dy);
                end;
            case 'gauss'
                for a = 1:4
                    Fd(a,ii,jj) = 1/4/pi*triplequad(@(x,y,z) ...
                        funG(x,y,z,xt,yt, a, dx, dy, h), ...
                        0, dx, 0, dy, -5*h, 5*h);
                end;
        end;
        counter_1 = counter_1 + 1;
        fprintf('.');
        if not(mod(counter_1,32))
            fprintf('\n');
        end;
    end;
end;
fprintf('\n');
%%% End of function

function Fd=FdMatrix_rect (nx,ny, dx, dy, h, zprofile)
n = nx*ny;
m = (nx-1)*(ny-1);
Fd = zeros(4,n,m);
F_temp = zeros(4,2*(nx-1),2*(ny-1));
counter_1 = 0;
for xp=1:2*(nx-1)
    for yp=1:2*(ny-1)
        xt = (xp - nx + 1)*dx; % For rectangular grid of electrodes
        yt = (yp - ny + 1)*dy;
        switch zprofile
            case 'step'
                for a = 1:4
                    F_temp(a,xp,yp) = 1/2/pi*dblquad(@(x,y) ...
                        fun(x,y,xt,yt, a, dx, dy, h), 0,dx,0,dy);
                end;
                case 'gauss'
                for a = 1:4
                    F_temp(a,xp,yp) = 1/4/pi*triplequad(@(x,y,z) ...
                        funG(x,y,z,xt,yt, a, dx, dy, h), ...
                        0, dx, 0, dy, -5*h, 5*h);
                end;
        end;
                
        counter_1 = counter_1 + 1;
        fprintf('.');
        if not(mod(counter_1,32))
            fprintf('\n');
        end;
    end;
end;
for jj=1:m
    for ii=1:n
        [xi,yi] = ind2sub([nx,ny],ii);
        [xj,yj] = ind2sub([nx-1,ny-1],jj);
        for a=1:4
            Fd(a,ii,jj) = F_temp(a,xi-xj+nx-1, yi-yj+ny-1);
        end;
    end;
end;
fprintf('\n');
%%% End of function

function f = fun(x, y, xt, yt, a, dx, dy, h)
dist = sqrt( (xt-x).^2+(yt-y)^2);
value = ones(1,length(dist));
switch a
    case 2
        value = x/dx;
    case 3
        value = y/dy;
    case 4
        value = x/dx.*y/dy;
end;
f = value.*asinh(h./dist);
% End of function

function f = funG(x, y, z, xt, yt, a, dx, dy, h)
dist = sqrt( (xt-x).^2+(yt-y)^2+z^2);
dist(dist<1e-8)=1e-8;
value = ones(1,length(dist));
switch a
    case 2
        value = x/dx;
    case 3
        value = y/dy;
    case 4
        value = x/dx.*y/dy;
end;
f = value./dist*exp(-z^2/(2*h^2));
% End of function

function E=EMatrix (nx,ny)
n = nx*ny;
m = (nx-1)*(ny-1);
E = zeros(4,m,n);
for ii=1:m
    [xk,yk] = ind2sub([nx-1,ny-1],ii);
    nr = sub2ind([nx ny],xk,yk);
    nrx = sub2ind([nx ny],xk+1,yk);
    nry = sub2ind([nx ny],xk,yk+1);
    nrxy = sub2ind([nx ny],xk+1,yk+1);

    E(1,ii,nr) = 1; 
    
    E(2,ii,nr) = -1;
    E(2,ii,nrx) = 1;
    
    E(3,ii,nr) = -1;
    E(3,ii,nry) = 1;
    
    E(4,ii,nr) = 1;
    E(4,ii,nrxy) = 1;
    E(4,ii,nrx) = -1;
    E(4,ii,nry) = -1;
end;
% End of function
