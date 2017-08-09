function [hfig] = plotEyeDistAnimation(DATA, HEADERS, figtitle, delaytime, accepted_validities, varargin)
%Function [hfig] = plotEyeDistAnimation(DATA, HEADERS, figtitle, delaytime, accepted_validities, varargin)
%
% Function plots the distance of the eyes of the participant read from the 
% DATA-matrix. figtitle specifies the name of the figure. Varargin may
% contain still images that are placed in the animation accordingly.
% Delaytime is to tune the loop slower if the computer is "too" fast.
% varargin format is what follows:
% {imagefile, coordinates, startrow, endrow}
% coordinates in the form [xstart xend ystart yend]

imageparameters = varargin;

% load imagefiles and do other tunings
for i=1:length(imageparameters)
    parameters = imageparameters{i};
    img1 = imread(parameters{1});
    
    % image goes up-and down, turn around
    for k=1:3
       img1(:, :, k) = flipud(img1(:, :, k));
    end
    img{i} = img1;
    
    % because of the structure of the loop later (begins from 2 rather
    % than 1), change the start-time 1 to 2 as it makes little difference
    % to the viewer and facilitates the loops
    if parameters{3} == 1
        parameters{3} = 2;
        imageparameters{i} = parameters;
    end
end


limits = [-100 100];

% search column values to use
xgl = colNum(HEADERS, 'XGazePosLeftEye');
ygl = colNum(HEADERS, 'YGazePosLeftEye');
xgr = colNum(HEADERS, 'XGazePosRightEye');
ygr = colNum(HEADERS, 'YGazePosRightEye');

xl = DATA{xgl};
yl = DATA{ygl};
xr = DATA{xgr};
yr = DATA{ygr};

% search distance values to use
rdc = colNum(HEADERS, 'DistanceRightEye');
ldc = colNum(HEADERS, 'DistanceLeftEye');

valr = DATA{colNum(HEADERS, 'ValidityRightEye')};
vall = DATA{colNum(HEADERS, 'ValidityLeftEye')};

valrp = (~ismember(valr, accepted_validities) * (-1000) ) + limits(1) - 0.01*limits(1);
vallp = (~ismember(vall, accepted_validities) * (-1000) ) + limits(1) - 0.04*limits(1);

mean_dist_l = mean(getColumn(DATA, ldc));
mean_dist_r = mean(getColumn(DATA, rdc));
common_mean = mean([mean_dist_l mean_dist_r]);

rowcount = rowCount(DATA);

starttime = getValue(DATA, 1, colNum(HEADERS, 'TETTime'));

axlimits = [0 1 0 1];

% create figure
scrsz = get(0,'ScreenSize');
hfig = figure('Position',[0.2*scrsz(3) 0.2*scrsz(4) scrsz(3)/2 scrsz(4)/1.5]);
%hfig = figure;
set(hfig, 'name', figtitle, 'numbertitle', 'off');



% initialize the upper axes for coordinates
a1 = subplot(2,2,1);
h1 = plot(xl(1), yl(1), 'o', xr(1), yr(1), 'x'); %, 'erasemode', 'background');
axis(axlimits);
%width = 0.40;
%set(gca, 'position', [(1-width)/2 0.5838 width 0.3412])
title('Gaze in the display');

% head position axis
a2 = subplot(2,2,2);
eyedist = 40; % mm
center = 400;
displaywidth = 510;

% find normal
normal_vector = cross([eyedist DATA{ldc}(1) - DATA{rdc}(1) 0], [0 0 1]);

h3 = plot(center, 15, 'o', [center-eyedist/2 center+eyedist/2], [DATA{ldc}(1) DATA{rdc}(1)], 'o', [center center-100*normal_vector(1)], ...
          [mean([DATA{ldc}(1) DATA{rdc}(1)]) mean([DATA{ldc}(1) DATA{rdc}(1)]) + 100*normal_vector(2)], [center-displaywidth/2 center+displaywidth/2], [15 15] );
set(gca, 'YDir','reverse');
axis([0 800 0 800]);
title(['Eye-positions and normal. Distance between eyes: ' num2str(eyedist) 'mm']);



% initialize the lower axes for time-view
a3 = subplot(2,2,3:4);
x = zeros(1, rowcount);
x(1) = 0;

% construct timevector
for i=2:rowcount
    x(i) = getValue(DATA, i, colNum(HEADERS, 'TETTime')) - starttime;
end

%x = 1:rowcount;
h2 = plot([0 0], [limits(1) limits(2)], 'r', x, DATA{ldc} - common_mean, x, DATA{rdc} - common_mean, x, valrp, '.', x, vallp, '.');
axis([min(x) max(x) limits(1) limits(2)]);

%xlabel(['Time: ' '0']);
title('Eye distances');
xlabel('Time from start (ms)');
ylabel('Normalized distance, mm');

original_ticks = get(gca,'YTickLabel');

for i=1:length(original_ticks)
    if (str2num(original_ticks(i,:)) == 0)
        new_ticks{i} = num2str(common_mean);
    else
        new_ticks{i} = original_ticks(i,:);
    end
end

set(gca, 'YTickLabel', new_ticks);

set(gcf, 'currentaxes', a1);

for i=2:rowcount

    % go through the images and plot those that have been selected
    for j=1:length(imageparameters)
        parameters = imageparameters{j};
        if parameters{3} == i
            hold on;            
            coord = parameters{2};
            siz = size(img{j});

            imx = [0:1:siz(2)]./siz(2).*(coord(2)-coord(1))+ coord(1);%-1/2*siz(2)/xres+0.5;
            imy = [0:1:siz(1)]./siz(1).*(coord(4)-coord(3))+ coord(3);%-1/2*siz(1)/yres+0.5;
            imghandle(j) = image(imx, imy, img{j}, 'parent', a1);
            uistack(imghandle(j), 'bottom');
            hold off;
        end
        
        if parameters{4} == i
            delete(imghandle(j));
        end
    end
    
    set(h1(1), 'Xdata', xl(1:i), 'Ydata', yl(1:i));
    set(h1(2), 'Xdata', xr(1:i), 'Ydata', yr(1:i));
    
    
    normal_vector = cross([eyedist DATA{ldc}(i) - DATA{rdc}(i) 0], [0 0 1]);
    set(h3(2), 'Ydata', [DATA{ldc}(i) DATA{rdc}(i)]);
    set(h3(3), 'Xdata', [center center-100*normal_vector(1)]);
    set(h3(3), 'Ydata', [mean([DATA{ldc}(i) DATA{rdc}(i)]) mean([DATA{ldc}(i) DATA{rdc}(i)]) + 100*normal_vector(2)]);
    
    
    set(h2(1), 'Xdata', [x(i) x(i)], 'Ydata', [limits(1) limits(2)]);
    %xlabel(a2, ['Time: ' num2str(round(getValue(DATA, i, colNum(HEADERS, 'TETTime')) - starttime))]);
    
    drawnow;
    
    pause(delaytime);
end