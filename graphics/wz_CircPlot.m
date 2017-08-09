function CIRC = wz_CircPlot(ang, w, nbins, meth, nopts, maxRad, jitt)
% wz_CircPlot - plot polar plot as (weighted) histogram
%
%  
% DESCRIPTION 
%
%  
% SYNTAX 
%   CIRC = wz_CircPlot(ang, w, nbins, nopts)
%
%   Input:
%           -ang    vector of angles associated with data points (in degree) 
%       
%           -w      weight for the data point (e.g. firing rate at a given direction)
%       
%           -nbins  number of bins for the polar histogram
%       
%           -meth   method to determine center location and error bars. per
%                   default, a normal approach is used, i.e. mean +/- ste
%                   options are: 'mean' (default), 'median', 'poisson'
%       
%           -nopts  do not plot individual data points, just use 
%       
%           -maxRad define maximal scale range
%       
%           -jitt   jitter data points
%       
%       
%   
%   Output:
%       
%
% REFERENCES 
%
%
% ......................................................................... 
% wolf zinke, wolfzinke@gmail.com 
%
% $Created : 06-Oct-2014 by wolf zinke
% $Modified: 

% ____________________________________________________________________________%
%% check arguments and define default values 

if(~exist('w','var') || isempty(w))
    w = ones(size(ang));
elseif(length(w) == 1)
    w = repmat(w,1,length(ang));
elseif(length(ang) ~= length(w))
    error('Both vectors must have the same length!');
end

neg = ang < 0;
if(any(neg) == 1)
    ang(neg) = ang(neg) + 360;
end

if(~exist('meth','var') || isempty(meth))
    meth = 'mean';
end
if(~exist('nbins','var') || isempty(nbins))
    nbins = 8;
end

if(~exist('nopts','var') || isempty(nopts))
    nopts = 0;
end

if(~exist('maxRad','var'))
    maxRad = [];
end

if(~exist('jitt','var') || isempty(jitt))
    jitt  = 0;
end

% remove nan values
ang = rowcheck(ang);
w   = rowcheck(w);

p = ~isfinite(ang) | ~isfinite(w);
ang(p) = [];
w(p)   = [];

% ____________________________________________________________________________%
%% process data
radval = deg2rad(ang);

angcntr = [0: 360/nbins : 360-360/nbins];
angbins = angcntr -(180/nbins);
angbins(end+1) = 360 +(180/nbins);

seccnts = nan(size(angcntr));
if(length(unique(w)) > 1)
    secvals = nan(size(angcntr));
    secCIl  = nan(size(angcntr));
    secCIu  = nan(size(angcntr));
end


for(i=1:length(angbins)-1)
    if(angbins(i)<0)
        pp = ang >= angbins(i)+360 | ang < angbins(i+1);
    else
        pp = ang >= angbins(i) & ang < angbins(i+1);
    end
    
    seccnts(i) = sum(pp);
    if(sum(pp) > 0)
        switch meth
            case 'mean'
                secvals(i) = mean(w(pp));
                secste     = std(w(pp)) / sqrt(seccnts(i));
                secCIl(i)  = secvals(i) - secste;
                secCIu(i)  = secvals(i) + secste;
                
            case 'median'
                secvals(i) = median(w(pp));
                
 %               secste = 1.57*diff(prctile(w(pp),[25,75]))/sqrt(seccnts(i));
                secste = 1.4826*mad(w(pp),1) / sqrt(seccnts(i));
                secCIl(i)  = secvals(i) - secste;
                secCIu(i)  = secvals(i) + secste;

            case 'poisson'
                [secvals(i),lambdaci] = poissfit(w(pp));
                secvals(i) = poisstat(secvals(i));
                secCIl(i)  = lambdaci(1);
                secCIu(i)  = lambdaci(2);
                
            otherwise
                error('Method unknown!');
        end
    end
end

% remove missing data 
novals = seccnts == 0 | isnan(secvals);
seccnts(novals)  = [];
secvals(novals)  = [];
secCIl(novals)   = [];
secCIu(novals)   = [];
angcntr(novals)  = [];

if(sum(novals) > 0)
    warning('Not enough data for each angle bin. Inacurate results!');
end

radcntr = deg2rad(angcntr);

angstep = unique(diff(radcntr));

if(length(angstep > 1))
    angstep = [];
end

% get some summary statistics
    CIRC.angles    = ang;
    CIRC.weights   = w;
    CIRC.angcntr   = angcntr;
    CIRC.method    = meth;
    CIRC.secvals   = secvals;
    CIRC.secCIl    = secCIl;
    CIRC.secCIu    = secCIu;
    
    CIRC.Nbins     = length(secvals); 
    CIRC.EquiDist  = length(unique(diff(angcntr))) == 1; % are data points equally spaced?
    
    % summary statistics based on raw data
    if(length(radval) > 4)
        CIRC.MU_raw      = circ_mean(radval(:), w(:)); % mean angle calculate over full raw data
        CIRC.std_raw     = circ_std(radval(:),  w(:));  % circular standard deviation of raw data
        CIRC.AngMean_raw = rad2deg(CIRC.MU_raw);
        CIRC.AngStd_raw  = rad2deg(CIRC.std_raw);
        CIRC.R_raw       = circ_r(radval(:),    w(:));    % mean vector length
        CIRC.Conf_raw    = circ_confmean(radval(:), 0.05, w(:));  % 95% confidence interval
%       CIRC.Otest_raw   = circ_otest(radval(:), 1, w(:));
        [CIRC.Vtest_raw, CIRC.VtestV_raw] = circ_vtest(radval(:), CIRC.MU_raw, w(:));
        [CIRC.Rtest_raw, CIRC.RtestZ_raw] = circ_rtest(radval(:), w(:));
        
        if(CIRC.AngMean_raw < 0)
            CIRC.AngMean_raw = 360 + CIRC.AngMean_raw;
        end
    else
        CIRC.MU_raw      = NaN;
        CIRC.std_raw     = NaN;
        CIRC.AngMean_raw = NaN;
        CIRC.AngStd_raw  = NaN;
        CIRC.R_raw       = NaN;
        CIRC.Conf_raw    = NaN;
        CIRC.Vtest_raw   = NaN;
        CIRC.VtestV_raw  = NaN;
%        CIRC.Otest_raw  = NaN;
        CIRC.Rtest_raw   = NaN;
        IRC.RtestZ_raw   = NaN;
    end
    
    % summary statistics based on binned dat
    if(length(secvals) > 2)
        CIRC.MU_bin      = circ_mean(radcntr(:), secvals(:));  % mean angle calculate over bins
        CIRC.std_bin     = circ_std(radcntr(:),  secvals(:));   % circular standard deviation of binned data
        CIRC.AngMean_bin = rad2deg(CIRC.MU_bin);
        CIRC.AngStd_bin  = rad2deg(CIRC.std_bin);
        CIRC.R_bin       = circ_r(radcntr(:),    secvals(:), angstep); % mean vector length
        CIRC.Conf_bin    = circ_confmean(radcntr(:), 0.05, secvals(:));  % lower bound of 95% confidence interval
%       CIRC.Otest_bin   = circ_otest(radcntr(:), intStep, secvals(:));
        [CIRC.Vtest_bin, CIRC.VtestV_bin] = circ_vtest(radcntr(:), CIRC.MU_bin, secvals(:));
        [CIRC.Rtest_bin, CIRC.RtestZ_bin] = circ_rtest(radcntr(:), secvals(:));

        CIRC.tune      = circ_r(radcntr(:), secvals(:)/max(secvals), angstep); % normalized vector length
        
        if(CIRC.AngMean_bin < 0)
            CIRC.AngMean_bin = 360 + CIRC.AngMean_bin;
        end
    else
        CIRC.MU_bin      = NaN;
        CIRC.std_bin     = NaN;
        CIRC.AngMean_bin = NaN;
        CIRC.AngStd_bin  = NaN;
        CIRC.R_bin       = NaN;
        CIRC.Conf_bin    = NaN;
        CIRC.Vtest_bin   = NaN;
        CIRC.VtestV_bin  = NaN;
%        CIRC.Otest_bin  = NaN;
        CIRC.Rtest_bin   = NaN;
        CIRC.RtestZ_bin  = NaN;

        CIRC.tune      = NaN;
    end
    
% ____________________________________________________________________________%
%% plot data
% plot data points

% set the axis limits
if(~isempty(maxRad))
    pw = w(w<maxRad);
else
    pw = w;
    if(nopts == 1)
        if(length(unique(pw)) > 1)
            maxRad = max(secCIu);
        else
            maxRad = max(secvals);
        end
    else
        maxRad = prctile(pw,90);
    end
    if(isempty(maxRad))
        maxRad = 1; % just set some value to create the plot
    end
end

% initialize plotting area
t = 0 : .01 : 2 * pi;
P = polar(t, maxRad * ones(size(t)));
set(P, 'Visible', 'off')
hold on

if(length(unique(w)) > 1 && nopts == 0)
    
    if(jitt == 1)
        plrad = radval + deg2rad((randi(50,size(radval))-25)./10);
    else
        plrad = radval;
    end
    polar(plrad, pw,'.k');
end

if(length(unique(pw)) > 1)
    % secplot(radcntr,secvals,'color','r','LineWidth',2);
    wz_polar(angcntr,  secvals, 'b', [], 'LineWidth',2);
%     h1 = polar(radcntr,  secvals, '-b');
%     set(h1,'LineWidth',2);
    h = polar(radcntr, secvals, 'ob');
    set(h, 'MarkerFaceColor', 'b')
    
%     pH = polar(radcntr,secvals,'.b');
    for(i=1:length(radcntr))
        cH = polar([radcntr(i), radcntr(i)], [secCIl(i), secCIu(i)]);
        set(cH,'color','b','LineWidth',2);
        
        text(1.4*maxRad*cos(radcntr(i)),1.4*maxRad*sin(radcntr(i)),int2str(seccnts(i)),...
            'HorizontalAlignment','center', 'VerticalAlignment','middle','color','b');
    end
    rH = polar([CIRC.MU_bin,CIRC.MU_bin], [0 CIRC.R_bin*max(secvals)]);
    set(rH,'color','r','LineWidth',2);
else
    secplot(deg2rad(angcntr),seccnts,'LineWidth',2);  
end
    
hold off    

% ____________________________________________________________________________%
%% helper function
function vec = rowcheck(vec)
    if size(vec, 1) ~= 1
        vec = vec';
    end

function varargout = secplot(th, r, varargin)
% ax = secplot(th, rho) Plot th and r as a series of triangular wedges
% 360/size(th) degrees wide, with a radius corresponding to r.
%
% Return a handle to the polar plot.
%
% Additional arguments will be passed to polar, which handles the plotting.
%
% Example:
%       th = 0:pi/4:2*pi;
%       secplot(th, sin(th/10));
%
% Matt Foster <ee1mpf@bath.ac.uk>
%
% $Id$

    % Sanity checking:
    if ~isvector(th) || ~isvector(r)
        warning('Secplot:ArgWarning', 'Inputs th and r should be vectors');
    end

    % These _must_ be rows -- make sure they are.
    th = rowcheck(th);
    r = rowcheck(r);

    % Find triangle spacing
    space = diff(th)./2;

    if any(diff(space) > 1e-6)
        space = mean(space);
        warning('SecPlot:ArgError', ...
        'input angles not evenly spaced, using mean value: %2.4g', ...
        space);
    else
        space = space(1);
    end

    % Calculate lower and upper bounds of triangles
    lw = th - space;
    up = th + space;

    % Repeat angles (each radial side has two of the same angle)
    angs = [lw; lw; up; up];
    angs = angs(:);

    % handle the radii (the non-radial sides consist of two points at each
    % radius)
    rads = [zeros(size(r)); r; r; zeros(size(r))];
    rads = rads(:);

    % Plot:
    ax = polar(angs, rads);
    set(ax, varargin{:});

    if nargout ~= 0
        varargout{:} = ax;
    end


