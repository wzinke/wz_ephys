function [h, thv, lpt] = wz_polar(ang, w, colstyl, ip, varargin)
% correct the linear interpolation in the polar plot function.
% Matlab is drawing polar plots the wrong way, it interpolates two
% points in the polar coordinate system with a straight line. This implies
% the wrong values, because the points of a line should be polar-transformed
% as well, resulting in a curve.
%
% wolf zinke, June 2014
%
% ToDo:  - add polygone functionality
%        - allow data points to be plotted 

if(~exist('ip','var') || isempty(ip)) % angle resolution for interpolation
    ip = 1;
end

if(~exist('colstyl','var') || isempty(colstyl))
    colstyl = '-b';
end

if(length(w) == 1)
    w = repmat(w, 1, length(ang));
elseif(length(ang) ~= length(w))
    error('Both vectors must have the same length!');
end

[~, idx] = sort(ang);
ang = ang(idx);
w   =   w(idx);

lpt = [];
apt = [];

% interpolate data between angle positions
for(i=1:length(ang)-1)
    angstp = diff(ang([i,i+1]))/ip;
    apt    = [apt, ang(i):ip:ang(i+1)-ip];

    if(diff(w([i,i+1])) == 0)
%     stp    = diff(w([i,i+1]))/angstp;
%     if(stp == 0)
        lpt = [lpt, repmat(w(i),1,angstp)];
    else
        % lpt = [lpt, w(i):stp:w(i+1)-stp];
        lpt = [lpt, linspace( w(i),w(i+1),angstp)];
    end
end

% connect now first with last point
angstp = (360-ang(end)+ang(1))/ip;
apt    = [apt, mod(ang(end):ip:(360+ang(1)), 360)];

if(diff(w([end,1])) == 0)
    lpt = [lpt, repmat(w(end),1,angstp+1)];
% stp    = diff(w([end,1]))/angstp;
% 
% if(stp == 0)
%     lpt = [lpt, repmat(w(end),1,ip)];
else
    %lpt = [lpt,w(end):stp:w(1)];
    lpt = [lpt, linspace( w(end),w(1),angstp+1)];
end

thv = deg2rad(apt); % ./ (180/pi);  % needs to be fixed and put into the loop

% if(length(thv) == length(lpt)-1)
%     lpt(end) == [];
% end

h = polar(thv, lpt, colstyl);

if(length(varargin) > 1)
    for(i=1:2:length(varargin))
        set(h,varargin{i},varargin{i+1});
    end
end


