function [intensity_normalized, xin, yin] = calculateIntensityMatrix(xcol, ycol, HM_xresolution, HM_yresolution)
%Function [intensity_normalized, xin, yin] = calculateIntensityMatrix(xcol, ycol, HM_xresolution, HM_yresolution)
%
% Returns intensity matrix calculated according to the gaze point spatial
% distribution. So the points where gaze spent most of the time appear with
% higher intensities on the matrix. HM_xresolution and HM_yresolution are
% number parameters which specify the resolution of the heatmap. Each point
% is rounded to the following resolution "pixel". Coordinates for heatmap
% are returned in xin, yin. Coordinates are centered "half-step" from zero
% to one minus half step according to the step being 1/HMresolution.
% Intensity heatmap values are also normalized by dividing each point's
% intensity value by total number of samples.

% calculate and plot heatmap over the imagelayer
disp('Calculating gaze intensity matrix...');

% form intensity-matrix (+1 to add zero)
intensity = zeros(HM_yresolution, HM_xresolution);

% round datapoints to match a unique point in the matrix

% go through all the datapoints
dpoints = 0;
for i=1:length(xcol)

    rx = ceil(xcol(i)*HM_xresolution);
    ry = ceil(ycol(i)*HM_yresolution);

    % the edges of the screen
    if rx==0
        rx=1;
    elseif ry==0
        ry=1;
    end
    
    if (1 <= rx && rx <= HM_xresolution) && (1 <= ry && ry <= HM_yresolution)
        % increase the the intensity of the nearest point if inside screen area
        intensity(ry, rx) = intensity(ry, rx) + 1;
        dpoints = dpoints + 1;
    end
end

% normalize intensity matrix
%intensity_normalized = intensity./max(max(intensity));
intensity_normalized = intensity./(dpoints);

% Make a truecolor all-green image. (if you want one color-image)
%green = cat(3, zeros(size(intensity)), ones(size(intensity)), zeros(size(intensity)));

xi=0:HM_xresolution-1;
xi = xi+0.5;
xin = xi./HM_xresolution;
yi=0:HM_yresolution-1;
yi = yi+0.5;
yin = yi./HM_yresolution;

disp('Done.');