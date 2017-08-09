%   Copyright 2006-2010 Szymon Leski
%   s.leski@nencki.gov.pl

function csd_out = Get_csd(name, method, boundary, h)

% method = 'trad', 'step', 'lin', 'spline' (not-a-knot)
% boundary = 'no'='S', 'B', 'D', 'Vaknin'

datadir = 'data';
datasetdir = 'datasets';
if strcmp(boundary, 'S')
    boundary = 'no'; % so the file names are consistent
end; % -> boundary = 'no', 'B', 'D', 'Vaknin'

% Load grid info
potsname = fullfile(datasetdir, [name '-pots.mat']);
data = load(potsname, 'dx', 'dy', 'nx', 'ny');
VX = 1:0.05:data.nx;
VY = 1:0.05*(data.dx/data.dy):data.ny;

% Compose the csd file name
[dynum, dyden] = rat(data.dy/data.dx);
switch method
    case 'trad'
        hnum = 1;
        hden = 1;
        nx2 = data.nx;
        ny2 = data.ny;
    otherwise
        [hnum, hden] = rat(h/data.dx);
        switch boundary
            case 'no'
                nx2 = data.nx;
                ny2 = data.ny;
            case {'B', 'D'}
                nx2 = data.nx + 2;
                ny2 = data.ny + 2;
        end
end
Fname = [num2str(nx2) '-' num2str(ny2) '-' ...
    num2str(dynum) '-' num2str(dyden) '-' ...
    num2str(hnum) '-' num2str(hden)  '-' method]; % czasem nie zdef nx2 ny2
csdname = fullfile(datasetdir, [name '-' ...
    method '-' boundary '-' ...
    num2str(hnum) '-' num2str(hden)  '.mat']);
% Check if cached CSD is available
if exist(csdname, 'file')
    % If yes, load and return
    disp('  Loading CSD')
    csddata = load(csdname, 'csd2');
    csd_out = csddata.csd2;
else
    % If no, calculate, interpolate, save, and return
    % Load potentials
    data = load(potsname, 'dx', 'dy', 'nx', 'ny', 'pots');
    switch method
        case 'trad'
            switch boundary
                case 'no'
                    Finv = inittrad2d(data.nx, data.ny, data.dx, data.dy);
                case 'Vaknin'
                    Finv = inittrad2d(data.nx, data.ny, ...
                        data.dx, data.dy, 'Vaknin');
                otherwise
                    error(['Boundary ' boundary ' illegal for trad csd.']);
            end;
            disp('  Calculating CSD')
            csd = icsd2d(data.pots, Finv, 'trad');
            disp('  Interpolating CSD')
            csd2 = interp2d(csd, VX, VY, 'splinem'); 
        case {'step', 'lin', 'spline'}
            % Check if precalculated F matrix is available
            % If not, create it
            if ~exist(fullfile(datadir,[Fname '.mat']), 'file')
                disp('  Creating CSD matrices');
                % rescale dx to 1, ie. divide all lengths by dx
                switch method
                    case 'step'
                        initstep2d(Fname, nx2, ny2, 1, ...
                            data.dy/data.dx, h/data.dx);
                    case 'lin'
                        initlin2d(Fname, nx2, ny2, 1, ...
                            data.dy/data.dx, h/data.dx);
                    case 'spline'
                        initspline2d(Fname, nx2, ny2, 1, ...
                            data.dy/data.dx, h/data.dx);
                end
            end
            disp('  Loading CSD matrices');
            switch method
                case {'step', 'lin'}
                    Fdata = load(fullfile(datadir,[Fname '.mat']), 'F');
                    F = Fdata.F;
                case 'spline'
                    Fdata = load(fullfile(datadir,[Fname '.mat']), 'Fm');
                    F = Fdata.Fm;
            end
            % Calculate CSD
            disp('  Calculating CSD')
            csd = icsd2d(data.pots, F, boundary)/(data.dx^2);
            disp('  Interpolating CSD')
            switch method
                case {'step', 'lin'}
                    csd2 = interp2d(csd, VX, VY, ...
                        method, boundary); %#ok<NASGU>
                case 'spline' % = splinem
                    csd2 = interp2d(csd, VX, VY, ...
                        'splinem', boundary); %#ok<NASGU>
            end
        otherwise
            error(['Method ' method ' not recognized']);
    end
    csd_out = csd2;
    save(csdname, 'csd2', 'csd');
end
disp('  Done with csd.')
