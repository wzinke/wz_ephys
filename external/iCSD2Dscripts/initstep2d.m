function initstep2d(name,nx,ny,dx,dy,h,zprofile,dsp)

% INITSTEP2D Calculates F matrix for 2-D step ICSD method.
%
%   initstep2d(name, nx, ny, dx, dy, h)
%   initstep2d(name, nx, ny, dx, dy, h, dsp)
%   initstep2d(name, nx, ny, dx, dy, h, zprofile)
%   initstep2d(name, nx, ny, dx, dy, h, zprofile, dsp)
%
%   Calculates the matrix F which relates CSD to potentials.
%   F(k,l) is the potential generated at k-th node by 
%   unit CSD placed at l-th node of 2-D grid if nearest
%   neighbor interpolation of CSD between nodes is assumed.
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

parseinpargs

disp(' Generating F matrix for step distribution...');
if rectgrid
    F = FdMatrix_rect(nx, ny, dx, dy, h, zprofile); 
    % faster, rectangular only
    save(fname, 'F', 'name', 'dx', 'dy', 'h', ...
        'nx', 'ny', 'zprofile');
else
    F = FdMatrix(nx, ny, dx, dy, h, zprofile, dsp); 
    % slower, general geometry
    save(fname, 'F', 'name', 'dx', 'dy', 'h', ...
        'nx', 'ny', 'dsp', 'zprofile');
end;

function Fd=FdMatrix (nx,ny, dx, dy, h, zprofile, dsp)
n = nx*ny;
Fd = zeros(n);
counter_1 = 0;
for jj=1:n
    for ii=1:n
        [xi,yi] = ind2sub([nx,ny],ii);
        [xj,yj] = ind2sub([nx,ny],jj);
        xt = (xi - xj)*dx-dsp(xi,yi,1); % For non-rectangular 
        yt = (yi - yj)*dy-dsp(xi,yi,2); % grid of electrodes
        switch zprofile
            case 'step'
                Fd(ii,jj) = 1/2/pi*dblquad(@(x,y) ...
                    fun(x, y, xt, yt, h), -dx/2,dx/2,-dy/2,dy/2);
            case 'gauss'
                Fd(ii,jj) = 1/4/pi*triplequad(@(x,y,z) ...
                    funG(x, y, z, xt, yt, h), ...
                    -dx/2,dx/2,-dy/2,dy/2, -5*h, 5*h);
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
Fd = zeros(n);
F_temp = zeros(2*nx-1,2*ny-1);
counter_1 = 0;
for xp=1:2*nx-1
    for yp=1:2*ny-1
        xt = (xp - nx)*dx; % For rectangular grid of electrodes
        yt = (yp - ny)*dy;
        switch zprofile
            case 'step'
                F_temp(xp,yp) = 1/2/pi*dblquad(@(x,y) ...
                    fun(x, y, xt, yt, h), -dx/2,dx/2,-dy/2,dy/2);
            case 'gauss'
                F_temp(xp,yp) = 1/4/pi*triplequad(@(x,y,z) ...
                    funG(x, y, z, xt, yt, h), ...
                    -dx/2,dx/2,-dy/2,dy/2, -5*h, 5*h);
        end;
        counter_1 = counter_1 + 1;
        fprintf('.');
        if not(mod(counter_1,32))
            fprintf('\n');
        end;
    end;
end;
for jj=1:n
    for ii=1:n
        [xi,yi] = ind2sub([nx,ny],ii);
        [xj,yj] = ind2sub([nx,ny],jj);
        Fd(ii,jj) = F_temp(xi-xj+nx, yi-yj+ny);
    end;
end;
fprintf('\n');
%%% End of function

function f = fun(x, y, xt, yt, h)
dist = sqrt( (xt-x).^2+(yt-y)^2);
dist(dist<1e-8)=1e-8;
f = asinh(h./dist);
% End of function

function f = funG(x, y, z, xt, yt, h)
dist = sqrt( (xt-x).^2+(yt-y)^2+z^2);
dist(dist<1e-8)=1e-8;
f = 1./dist*(exp(-z^2/(2*h^2));
% End of function
