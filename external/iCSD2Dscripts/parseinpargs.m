% PARSEINPARGS Parse input arguments for CSD scripts. 

%   Copyright 2006-2010 Szymon Leski
%   s.leski@nencki.gov.pl

fname = fullfile(datafolder, [name '.mat']);
if exist(fname, 'file')
    error(['File ' fname ' already exists.']);
end;
if nargin==6            % Default - rectangular grid & step profile
    zprofile = 'step';
    rectgrid = 1;           
elseif nargin==7
    if strcmp(zprofile,'step')||strcmp(zprofile,'gauss')
        % Rect grid, profile specified
        rectgrid = 1;
    else % Non-rect grid, default (step) profile
        rectgrid = 0;
        dsp = zprofile;
        zprofile = 'step';
    end;
elseif nargin==8                % Non-rect grid, profile specified
    rectgrid = 0;           
else
    error('Wrong number of arguments')
end;
