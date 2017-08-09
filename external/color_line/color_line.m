function h = color_line(x,y,z,cmap)
%COLOR_LINE  Plot a vector with colors defined by a colormap.
%
% COLOR_LINE is a basic wrapper function which applies the matlab LINE
%            command, plotting a line for each {x, y} data pair with a
%            color set by the corresponding intensity value. The function
%            returns a column vector of handles to LINE objects.
%
%   Syntax:
%      handles = COLOR_LINE(X,Y,Z,CMAP)
%
%
%   Inputs:
%
%      X, Y       Spatial coordinate data
%      Z          Color data (intensity)
%      CMAP       Colormap, N x 3 array of RGB values
%
%
%   Example:
%               color_line( ...
%                          [-1:0.1:1]                   , ...
%                          [-1:0.1:1]                   , ...
%                          randn(length([-1:0.1:1]),1)  , ...
%                          'jet'                        )
%
% See also LINE, COLORMAP

if nargin < 4
    cmap = colormap('jet');
end

z_min = min(z);
z_max = max(z);

colormap(cmap);
cm = colormap;
cm_length = length(cm);

range = z_max - z_min;
cm_stepsize = range/cm_length;

h_list = zeros(length(x),1); %initialize list to store line object handles

hold on
for count=1:length(x)
    if ~isnan(z(count))
        cm_index = floor(...
            (z(count) - z_min) / ( cm_stepsize )...
            ) + 1;
        try %need this in case z_min is zero
            mc = cm(cm_index,:);
        catch
            mc = cm(cm_index-1,:);
        end
        h_list(count) = line(x(count),y(count)          ,...
            'LineWidth'             ,       1           ,...
            'Marker'                ,       'o'         ,...
            'MarkerEdgeColor'       ,       mc          ,...
            'MarkerFaceColor'       ,       mc          ,...
            'MarkerSize'            ,       10          ,...
            'color'                 ,       mc          );
    end
end
hold off

h = h_list;