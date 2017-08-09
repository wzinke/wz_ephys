function initspline2d(name,nx,ny,dx,dy,h,zprofile,dsp)

% INITSPLINE2D Calculates F matrix for 2-D spline ICSD method.
%
%   initspline2d(name, nx, ny, dx, dy, h)
%   initspline2d(name, nx, ny, dx, dy, h, dsp)
%   initspline2d(name, nx, ny, dx, dy, h, zprofile)
%   initspline2d(name, nx, ny, dx, dy, h, zprofile, dsp)
%
%   Calculates the matrix F which relates CSD to potentials.
%   F(k,l) is the potential generated at k-th node by 
%   unit CSD placed at l-th node of 2-D grid if spline
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
%   a .mat file in subdirectory data,
%   containing Fn, Fm (for natural and not-a-knot splines, respectively)
%   and some other variables

%   Copyright 2006-2010 Szymon Leski
%   s.leski@nencki.gov.pl

datafolder = 'data';

n = nx*ny;

parseinpargs

En = espmatrices2d(nx,ny,'n');
Em = espmatrices2d(nx,ny,'m');
disp(' Generating F matrix for spline distribution...');
if rectgrid
    F = FdSpMatrix_rect(nx, ny, dx, dy, h, zprofile); 
    % faster, rectangular only
else
    F = FdSpMatrix(nx, ny, dx, dy, h, zprofile, dsp); 
    % slower, general geometry
end;
Fn = zeros(n);
Fm = zeros(n);
for ii=1:4
    for jj=1:4
            Fn = Fn + squeeze(F(ii,jj,:,:))*squeeze(En(ii,jj,:,:));
            Fm = Fm + squeeze(F(ii,jj,:,:))*squeeze(Em(ii,jj,:,:));
    end;
end;
if rectgrid
    save(fname, 'Fn', 'Fm', 'name', 'dx', 'dy', 'h', ...
        'nx', 'ny', 'zprofile');
else
    save(fname, 'Fn', 'Fm', 'name', 'dx', 'dy', 'h', ...
        'nx', 'ny', 'dsp', 'zprofile');
end;

function Fd=FdSpMatrix (nx,ny, dx, dy, h, zprofile, dsp)
n = nx*ny;
m = (nx-1)*(ny-1);
Fd = zeros(4,4,n,m);
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
                    for b = 1:4
                        Fd(a,b,ii,jj) = 1/2/pi*dblquad(@(x,y) ...
                            funSp(x,y,xt,yt, a, b, dx, dy, h), 0,dx,0,dy);
                    end;
                end;
            case 'gauss'
                for a = 1:4
                    for b = 1:4
                        Fd(a,b,ii,jj) = 1/4/pi*triplequad(@(x,y,z) ...
                            funSpG(x,y,z,xt,yt, a, b, dx, dy, h), ...
                            0, dx, 0, dy, -5*h, 5*h);
                    end;
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

function Fd=FdSpMatrix_rect (nx,ny, dx, dy, h, zprofile)
n = nx*ny;
m = (nx-1)*(ny-1);
Fd = zeros(4,4,n,m);
F_temp = zeros(4,4,2*(nx-1),2*(ny-1));
counter_1 = 0;
for xp=1:2*(nx-1)
    for yp=1:2*(ny-1)
        xt = (xp - nx + 1)*dx; % For rectangular grid of electrodes
        yt = (yp - ny + 1)*dy;
        switch zprofile
            case 'step'
                for a = 1:4
                    for b = 1:4
                        F_temp(a,b,xp,yp) = 1/2/pi*dblquad(@(x,y) ...
                            funSp(x,y,xt,yt, a, b, dx, dy, h), 0,dx,0,dy);
                    end;
                end;
            case 'gauss'
                for a = 1:4
                    for b = 1:4
                        F_temp(a,b,xp,yp) = 1/4/pi*triplequad(@(x,y,z) ...
                            funSpG(x,y,z,xt,yt, a, b, dx, dy, h), ...
                            0, dx, 0, dy, -5*h, 5*h);
                    end;
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
for jj=1:m
    for ii=1:n
        [xi,yi] = ind2sub([nx,ny],ii);
        [xj,yj] = ind2sub([nx-1,ny-1],jj);
        for a=1:4
            for b=1:4
                Fd(a,b,ii,jj) = F_temp(a,b, xi-xj+nx-1, yi-yj+ny-1);
            end;
        end;
    end;
end;
%%% End of function

function f = funSp(x, y, xt, yt, a, b, dx, dy, h)
dist = sqrt( (xt-x).^2+(yt-y)^2);
value = ones(1,length(dist));

switch a
    case 1
        value = value.*(1-y/dy);
    case 2
        value = value.*y/dy;
    case 3
        value = 1/6*value.*((1-y/dy)^3-1+y/dy);%.*(dy^2);
    case 4
        value = 1/6*value.*((y/dy)^3-y/dy);%.*(dy^2);
end;
switch b
    case 1
        value = value.*(1-x/dx);
    case 2
        value = value.*x/dx;
    case 3
        value = 1/6*value.*((1-x/dx).^3-1+x/dx);%.*(dx^2);
    case 4
        value = 1/6*value.*((x/dx).^3-x/dx);%.*(dx^2);
end;
f = value.*asinh(h./dist);
% End of function

function f = funSpG(x, y, z, xt, yt, a, b, dx, dy, h)
dist = sqrt( (xt-x).^2+(yt-y)^2+z^2);
dist(dist<1e-8)=1e-8;
value = ones(1,length(dist));

switch a
    case 1
        value = value.*(1-y/dy);
    case 2
        value = value.*y/dy;
    case 3
        value = 1/6*value.*((1-y/dy)^3-1+y/dy);%.*(dy^2);
    case 4
        value = 1/6*value.*((y/dy)^3-y/dy);%.*(dy^2);
end;
switch b
    case 1
        value = value.*(1-x/dx);
    case 2
        value = value.*x/dx;
    case 3
        value = 1/6*value.*((1-x/dx).^3-1+x/dx);%.*(dx^2);
    case 4
        value = 1/6*value.*((x/dx).^3-x/dx);%.*(dx^2);
end;
f = value./dist*exp(-z^2/(2*h^2));
% End of function
