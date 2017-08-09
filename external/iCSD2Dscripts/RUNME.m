%   Copyright 2006-2010 Szymon Leski
%   s.leski@nencki.gov.pl

h = waitbar(0, 'Preparing data sets...');
load datasets/paper-data.mat
if ~exist(fullfile('datasets', 'subiculum-pots.mat'), 'file')
    Init_dataset('subiculum', subiculum.pots, ...
        subiculum.nx, subiculum.ny, ...
        subiculum.dx, subiculum.dy, subiculum.dt, subiculum.t0, ...
        subiculum.comment);
end;
waitbar(1/3, h, 'Preparing data sets...');
if ~exist(fullfile('datasets', 'model1-pots.mat'), 'file')
    Init_dataset('model1', model1.pots, model1.nx, model1.ny, ...
        model1.dx, model1.dy, model1.dt, model1.t0, model1.comment);
end;
waitbar(2/3, h, 'Preparing data sets...');
if ~exist(fullfile('datasets', 'model2-pots.mat'), 'file')
    Init_dataset('model2', model2.pots, model2.nx, model2.ny, ...
        model2.dx, model2.dy, model2.dt, model2.t0, model2.comment);
end;
close(h)

iCSD2Dtool
