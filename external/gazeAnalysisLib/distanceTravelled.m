function [distance_right, distance_left] = distanceTravelled(DATA, HEADERS, display_ratio, display_height)
%Function [distance_right, distance_left] = distanceTravelled(DATA, HEADERS, display_ratio, display_height)
%
% Calculates the distance moved by eyes during the DATA-section 
% and returns individual values for each eye. display ratio is a number
% that describes how many times wider display is than tall. E.G for 16:10
% display, use 1.6. display_height is the height of the display in SI unit,
% meters. E.G. for a 51cm tall display, use 0.51. EXPERIMENTAL.

% search column values to use
xgl = colNum(HEADERS, 'XGazePosLeftEye');
ygl = colNum(HEADERS, 'YGazePosLeftEye');
xgr = colNum(HEADERS, 'XGazePosRightEye');
ygr = colNum(HEADERS, 'YGazePosRightEye');

distance_right = 0;
distance_left = 0;

for i=2:rowCount(DATA)

    xldist = (DATA{xgl}(i-1) - DATA{xgl}(i)) * display_ratio * display_height;
    yldist = (DATA{ygl}(i-1) - DATA{ygl}(i)) * display_height;
    xrdist = (DATA{xgr}(i-1) - DATA{xgr}(i)) * display_ratio * display_height;
    yrdist = (DATA{ygr}(i-1) - DATA{ygr}(i)) * display_height;
    
    distance_right = distance_right + sqrt(xrdist^2+yrdist^2);
    distance_left = distance_left + sqrt(xldist^2+yldist^2);
end


%
% distance_right = xrdist*display_ratio*display_height + yrdist*display_height;
% distance_left = xldist*display_ratio*display_height + yldist*display_height;