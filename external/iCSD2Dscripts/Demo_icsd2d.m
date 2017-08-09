% DEMO_ICSD2D Demonstration of 2-D ICSD scripts

%   Copyright 2006-2010 Szymon Leski
%   s.leski@nencki.gov.pl

% Setup
datafolder = 'data';

% Generate F matrix for 8 x 16 electrodes array and linear interpolation
% The spacing of electrodes is 0.2 mm in x direction and 0.1 in y
% direction. Thickness of CSD distribution is 2h = 2x0.1 mm
nx = 8;
ny = 16;
dx = 0.2;
dy = 0.1;
h = 0.1;
if not(exist(fullfile(datafolder, 'example_lin.mat'), 'file'))
    initlin2d('example_lin', nx, ny, dx, dy, h);
end;
% For splines:
% initspline2d('example_spline', nx, ny, dx, dy, h);

% Load the calculated F matrix
load(fullfile(datafolder, 'example_lin'), 'F');
% For splines:
% load(fullfile(datafolder, 'example_spline'), 'Fn', 'Fm');


% Load example data (variable potentials).
% This is 1 x 8 x 16 array of potentials.
load(fullfile(datafolder, 'example.mat'));

% Calculate CSD
csd1 = icsd2d(potentials, F);
% For splines:
% csd1 = icsd2d(potentials, Fn); % natural splines
% csd1 = icsd2d(potentials, Fm); % not-a-knot splines (almost the same)

% Interpolate CSD to a dense grid
VX = 1:0.1:nx;
VY = 1:0.1:ny;
csd2 = interp2d(csd1, VX, VY, 'lin');
% For splines (natural)
% csd2 = interp2d(csd1, VX, VY, 'splinen'); % natural
% csd2 = interp2d(csd1, VX, VY, 'splinem'); % not-a-knot

% Plot the reconstructed CSD
clims = [-1 1]*max(abs(csd2(:)));
imagesc( (VX-1)*dx, (VY-1)*dy, squeeze(csd2(1,:,:))', clims);
axis image
colormap(flipud(jet))
