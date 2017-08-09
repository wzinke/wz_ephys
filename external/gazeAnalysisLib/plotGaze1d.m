function [hfig] = plotGaze1d(DATA, xgr, xgl, ygr, ygl)
%Function [hfig] = plotGaze1d(DATA, xgr, xgl, ygr, ygl)
%
% Plots the gaze in DATA one-dimensionally. A handle to figure is returned.
% xgr, ygr, etc. specify columns of gaze points for each eye to use.


disp('Plotting gaze (1d)...');

rowcount = rowCount(DATA);

hfig = figure;
plot([DATA{xgr} DATA{xgl} DATA{ygr} DATA{ygl}]);

axis([0 rowcount+1 0 1]);
%axis tight;

xlabel('Datapoint');
ylabel('Pixel coordinate');

disp('Done.');