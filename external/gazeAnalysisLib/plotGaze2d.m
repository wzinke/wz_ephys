function [hfig] = plotGaze2d(xg, yg)
% function [hfig] = plotGaze2d(xcol, ycol)
%
% Plots the gaze in the DATA-array 2-dimensionally. The figure handle hfig
% is returned. xg, yg are vectors of datapoints to plot.

disp('Plotting gaze (2d)...');

hfig = figure;
%plot(DATA{xgl}, DATA{ygl}, 'o', DATA{xgr}, DATA{ygr}, 'x');
plot(xg, yg, '.', 'markersize', 1);
set(gca,'YDir','reverse');
axis([0 1 0 1]);

disp('Done.');