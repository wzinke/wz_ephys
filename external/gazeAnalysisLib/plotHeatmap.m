function [hfig] = plotHeatmap(intensity_matrix, xhm, yhm, imageparameters, transparencycoefficient, transparencymode, clim, cmap, identifier)
%Function [hfig] = plotHeatmap(intensity_matrix, xhm, yhm, imageparameters, transparencycoefficient, transparencymode, clim, cmap, identifier)
%
% Plots the gaze in the DATA-array 2-dimensionally. The figure handle hfig
% is returned. intensity_matrix is the matrix containing intensity values,
% xhm and yhm the coordinates, imageparameters contains info about the
% images to display on the canvas, clim is color limits ([cmin cmax]),
% identifier the text below heatmap, cmap colormap, transparencymode (1,2)
% where 1=all that is zero is transparent and everything else not, 2=scaled
% transparency. Transparencycoefficient is a coefficient to scale
% transparency-level. For more detailed description of imageparameters,
% please see help plotGazeAnimation.

axlimits = [0 1 0 1];

% Load images
for i=1:length(imageparameters)
    parameters = imageparameters{i};
    img{i} = imread(parameters{1});
end

disp('Plotting gaze heatmap with images...');
hfig = figure;
set(hfig, 'numbertitle', 'off', 'name', identifier);
set(gca, 'YDir', 'reverse');
axis(axlimits);


% go through the images and plot those that have been selected
for i=1:length(imageparameters)
    parameters = imageparameters{i};

    hold on;
    coord = parameters{2};
    siz = size(img{i});

    imxn = [0:1:siz(2)]./siz(2).*(coord(2)-coord(1)) + coord(1);%-1/2*siz(2)/xres+0.5;
    imyn = [0:1:siz(1)]./siz(1).*(coord(4)-coord(3)) + coord(3);%-1/2*siz(1)/yres+0.5;
 
    imghandle(i) = image(imxn, imyn, img{i});
    set(imghandle(i),'CDataMapping','scaled');%otherwise,it will effect the colorbar
    uistack(imghandle(i), 'bottom');
    hold off;
end

% calculate and plot heatmap over the imagelayer
disp('Drawing gaze heatmap...');

hold on;
colormap(cmap)
%colormap('hot');
%colormap(flipud(colormap));
imagetoplot = intensity_matrix;
%imagetoplot = green;
himg = imagesc(xhm, yhm, imagetoplot);% *10000);
hold off;
set(gca, 'Clim', clim);
colorbar;


if transparencymode == 1 % all zero transparent, nothing else
    intensity_alpha = intensity_matrix;
    intensity_alpha(intensity_alpha ~= 0) = 1;
elseif transparencymode == 2 % scaled transparency
    intensity_alpha = intensity_matrix * transparencycoefficient;
    intensity_alpha(intensity_alpha > 1) = 1;
end

set(himg, 'alphadata', intensity_alpha); %,0.4);

disp('Done.');