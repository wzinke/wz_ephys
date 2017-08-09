function [distvector] = distanceBetwEyesInTheScreen(DATA, xgl, xgr, ygl, ygr)
%Function [distvector] = distanceBetwEyesInTheScreen(DATA, xgl, xgr, ygl, ygr)
%
% Returns the vector of distances between eyes in the screen. xgl, xgr, ygl,
% ygr are the corresponding numbers of the columns for xdata and ydata for
% both eyes.

disp('Calculating distance between eyes in the screen...');

rowcount = rowCount(DATA);

% calculate distance between x and y coordinate in 
for i=1:rowcount
    xdiff = abs(DATA{xgl}(i) - DATA{xgr}(i));
    ydiff = abs(DATA{ygl}(i) - DATA{ygr}(i));
    distvector(i) = sqrt(xdiff^2+ydiff^2);
end

disp('Done.');