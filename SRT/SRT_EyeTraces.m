function Gaze = SRT_EyeTraces(EyeX, EyeY, xvec, StimOn, GoCue, do_plot)
% Quick and dirty, data driven saccade detection for trial based data.
%
% DESCRIPTION
%
%
% SYNTAX
%
%   Gaze = SRT_EyeTraces(EyeX, EyeY, xvec, StimOn, GoCue, do_plot)
%
%   Input:
%
% .........................................................................
% wolf zinke, wolfzinke@gmail.com
%
% $Created : 19-Jun-2015 by wolf zinke
%


% ____________________________________________________________________________ %
%% parameter for plotting and saccade estimation
fixwin = 100; % time window before stimulus onset to determine fixation parameter
pltwin = [-1 * fixwin, fixwin];

mindur = 10;  % minimum duration above threshold to be considered a saccade
madthr = 6;   % how many multiples of the median absolute deviation define the threshold
percCI = [5 95]; % used to define the fixation window, i.e. the disperseness of the fixation

% ____________________________________________________________________________ %
%% check input arguments
if(~exist('StimOn','var') || isempty(StimOn))
    StimOn = 0; 
end

if(~exist('GoCue','var') || isempty(GoCue))
    GoCue = StimOn; % number of data points
end

if(~exist('do_plot','var') || isempty(do_plot))
    if(nargout >= 1)
        do_plot = 0;
    else
        do_plot = 1;
    end
end

if(~exist('xvec','var') || isempty(xvec))
    if(isvector(EyeX))
        xvec = 1:length(EyeX);
    else
        xvec = 1:size(EyeX,2);
    end
end

if(any(size(EyeX) ~= size(EyeY)))
    error('X and Y traces have different dimensions!');
end

% ____________________________________________________________________________ %
%% smooth data

smX = SRT_SmoothTrace(EyeX);
smY = SRT_SmoothTrace(EyeY);

% get the distance relative to previous eye position
Dist = getDist(EyeX, EyeY);
smD  = getDist(smX, smY);

% smD2 = SRT_SmoothTrace(Dist);

% ____________________________________________________________________________ %
%% determine relevant events

% determine fixation median coordinate as reference

% get pre-stimulus fixation
fpos = xvec >= StimOn - fixwin & xvec < StimOn;

fixX = nanmedian(smX(fpos));
fixY = nanmedian(smY(fpos));

madX = madthr*mad(smX(fpos), 1);
madY = madthr*mad(smY(fpos), 1);
madD = madthr*mad(smD(fpos), 1);

CI_X = prctile(smX(fpos), percCI); 
CI_Y = prctile(smY(fpos), percCI); 

noFix = smD > madD & xvec > StimOn;
saccstart = xvec(strfind(noFix, ones(1,mindur)));

if(~isempty(saccstart))
    saccstart = xvec(find(diff(smD) >= 0 & smD(2:end) > madD & xvec(2:end) > saccstart(1), 1, 'first')-1);
    saccend   = xvec(find(diff(smD) <= 0 & smD(2:end) < madD & xvec(2:end) > saccstart,    1, 'first')+1);
    if(isempty(saccend))
        saccend = NaN;
    end
else
    saccstart = NaN;
    saccend   = NaN;
end

if(isfinite(saccstart) && isfinite(saccend))
    ppp = xvec > saccstart & xvec < saccend;
    peak_val = max(smD(ppp));
    if(~isempty(peak_val))
        peak_time = xvec(find(smD == peak_val & ppp == 1, 1, 'first'));
    else
        peak_time = NaN;
    end
else
    peak_time = NaN;
%     peak_val  = NaN;    
end

% get pre-saccade fixation
ppos = xvec >= saccstart - fixwin & xvec < saccstart;

preX = nanmedian(smX(ppos));
preY = nanmedian(smY(ppos));

madXp = madthr*mad(smX(ppos), 1);
madYp = madthr*mad(smY(ppos), 1);
madDp = madthr*mad(smD(ppos), 1);

CIp_X = prctile(smX(ppos), percCI); 
CIp_Y = prctile(smY(ppos), percCI); 

% get post-saccade fixation
spos = xvec > saccend  & xvec < saccend+fixwin;

postX = nanmedian(smX(spos));
postY = nanmedian(smY(spos));

madXs = madthr*mad(smX(spos), 1);
madYs = madthr*mad(smY(spos), 1);
madDs = madthr*mad(smD(spos), 1);

CIs_X = prctile(smX(spos), percCI); 
CIs_Y = prctile(smY(spos), percCI); 

% ____________________________________________________________________________ %
%% define output struct
if(nargout >= 1)
    Gaze.EyeX = EyeX;
    Gaze.EyeY = EyeY;
    Gaze.smX  = smX;
    Gaze.smY  = smY;
    Gaze.Dist = Dist;
    Gaze.smD  = smD;
    
    Gaze.StimOn = StimOn;
    Gaze.GoCue  = GoCue;

    Gaze.SaccStart = saccstart;
    Gaze.SaccEnd   = saccend;
    Gaze.SaccDur   = saccend - saccstart;
    [Gaze.SaccAng,   Gaze.SaccAmp] = cart2pol(postX-preX,postY-preY);
    Gaze.SaccSpeed = Gaze.SaccAmp / Gaze.SaccDur;   
    
    saccpos = find(xvec >= saccstart & xvec <= saccend);
    if(~isempty(saccpos))
        Gaze.SaccDist   = sum(smD(saccpos));

        Gaze.Startpoint_X = smX(saccpos(1));
        Gaze.Startpoint_Y = smY(saccpos(1));

        Gaze.Endpoint_X = smX(saccpos(end));
        Gaze.Endpoint_Y = smY(saccpos(end));
    else
        Gaze.SaccDist     = NaN;
        Gaze.Startpoint_X = NaN;
        Gaze.Startpoint_Y = NaN;

        Gaze.Endpoint_X = NaN;
        Gaze.Endpoint_Y = NaN;
    end
    
    Gaze.fix_X = preX;
    Gaze.fix_Y = preY;
    Gaze.fix_madX = madX/madthr;
    Gaze.fix_madY = madY/madthr;
    Gaze.fix_madD = madD/madthr;
    Gaze.fix_CI_X = CI_X; 
    Gaze.fix_CI_Y = CI_Y;
    Gaze.fix_CI_Xrng = diff(CI_X);
    Gaze.fix_CI_Yrng = diff(CI_Y);

    Gaze.pre_X = fixX;
    Gaze.pre_Y = fixY;
    Gaze.pre_madX = madXp/madthr;
    Gaze.pre_madY = madYp/madthr;
    Gaze.pre_madD = madDp/madthr;
    Gaze.pre_CI_X = CIp_X; 
    Gaze.pre_CI_Y = CIp_Y;
    Gaze.pre_CI_Xrng = diff(CIp_X);
    Gaze.pre_CI_Yrng = diff(CIp_Y);

    Gaze.post_X = postX;
    Gaze.post_Y = postY;
    Gaze.post_madX = madXs/madthr;
    Gaze.post_madY = madYs/madthr;
    Gaze.post_madD = madDs/madthr;    
    Gaze.post_CI_X = CIs_X; 
    Gaze.post_CI_Y = CIs_Y;
    Gaze.post_CI_Xrng = diff(CI_X);
    Gaze.post_CI_Yrng = diff(CI_Y);
    
    if(saccstart <= GoCue)
        Gaze.Resp = 'early';
    elseif(isnan(saccstart))
        Gaze.Resp = 'miss';   
    else
        Gaze.Resp = 'decide';   
    end
end

% ____________________________________________________________________________ %
%% plot data
if(do_plot == 1)
    
    plpp = [ find(xvec > StimOn+pltwin(1), 1, 'first') : find(xvec > saccend+pltwin(2), 1, 'first') ];
    % plpp = [ find(xvec >= StimOn+pltwin(1), 1, 'first') : find(xvec >= StimOn & isfinite(EyeX) & isfinite(EyeY), 1, 'last') ];
    
    figure('Position', [0 0 1280 800], 'Renderer', 'Painters', 'Color', [1 1 1]);
    
    Yrng1 = minmax([EyeX(plpp);smX(plpp)]);
    Yrng2 = minmax([EyeY(plpp); smY(plpp)]);
    
    plyrng = 1.05 * max(abs([diff(Yrng1),diff(Yrng2)]));
    
    Xrng = [StimOn + pltwin(1) ; saccend + pltwin(2)];
    
    %%%%
    subplot('Position', [0.035 0.695 0.44 0.265]);
    hold on;
    set(gca,'TickDir','out','fontsize',10, 'XTick', []);
    plot(xvec(plpp), EyeX(plpp),'color',[0.4 0.4 0.4],'LineWidth',0.8);
    plot(xvec(plpp), smX(plpp),'b','LineWidth',2.5);
    
    Yrng = [-1, 1] * (plyrng/2) + mean(minmax(EyeX(plpp)));
    ylim(Yrng);
    
    xlim(Xrng);
    
    vline(StimOn,'color','g');
    vline(saccstart,'color','m', 'LineStyle','--');
    vline(saccend,'color','m', 'LineStyle','--');
    vline(peak_time, 'color', 'c', 'LineStyle',':');
    
    hline(fixX, 'color', 'g', 'LineStyle','-');
    hline(fixX+madX, 'color', 'g', 'LineStyle',':');
    hline(fixX-madX, 'color', 'g', 'LineStyle',':');

    hline(postX, 'color', 'r', 'LineStyle','-');
    hline(postX+madXs, 'color', 'r', 'LineStyle',':');
    hline(postX-madXs, 'color', 'r', 'LineStyle',':');
    
    %%%%
    subplot('Position', [0.035 0.385 0.44 0.265]);
    hold on;
    set(gca,'TickDir','out','fontsize',10, 'XTick', []);
    plot(xvec(plpp), smY(plpp),'b','LineWidth',2.5);
    plot(xvec(plpp), EyeY(plpp),'color',[0.4 0.4 0.4],'LineWidth',0.8);
    
    Yrng = [-1, 1] * (plyrng/2) + mean(minmax(EyeY(plpp)));
    ylim(Yrng);
    xlim(Xrng);
    
    vline(StimOn,    'color','g');
    vline(saccstart, 'color','m', 'LineStyle','--');
    vline(saccend,   'color','m', 'LineStyle','--');
    vline(peak_time, 'color', 'c', 'LineStyle',':');

    hline(fixY,      'color', 'g', 'LineStyle','-');
    hline(fixY+madY, 'color', 'g', 'LineStyle',':');
    hline(fixY-madY, 'color', 'g', 'LineStyle',':');

    hline(postY,       'color', 'r', 'LineStyle','-');
    hline(postY+madYs, 'color', 'r', 'LineStyle',':');
    hline(postY-madYs, 'color', 'r', 'LineStyle',':');


    %%%%
    subplot('Position', [0.035 0.065 0.44 0.265]);
    hold on;
    set(gca,'TickDir','out','fontsize',10);
    plot(xvec(plpp), smD(plpp),'b','LineWidth',2.5);
    plot(xvec, Dist,'color',[0.4 0.4 0.4],'LineWidth',0.8);
%    plot(xvec(plpp), smD2(plpp),'r','LineWidth',2.5);
    
    xlim(Xrng);
    
    vline(StimOn,    'color','g');
    vline(saccstart, 'color','m', 'LineStyle','--');
    vline(saccend,   'color','m', 'LineStyle','--');
    vline(peak_time, 'color', 'c', 'LineStyle',':');

    hline(madD,  'color', 'g', 'LineStyle',':');
    hline(madDs, 'color', 'r', 'LineStyle',':');
    
    xlabel('time [ms]');
    
    %%%%
    subplot('Position', [0.525 0.065 0.45 0.895]);
    hold on;
    
    plot(EyeX(plpp), EyeY(plpp), '.', 'color',[0.4 0.4 0.4])
    plot(smX(plpp), smY(plpp), '.b')
    
    ellipse(fixX,fixY,[diff(CI_X), diff(CI_Y)]./2,     'color','k', 'LineStyle',':', 'LineWidth',1.5);
    ellipse(preX,preY,[diff(CIp_X), diff(CIp_Y)]./2,   'color','g', 'LineStyle',':', 'LineWidth',1.5);
    ellipse(postX,postY,[diff(CIs_X), diff(CIs_Y)]./2, 'color','r', 'LineStyle',':', 'LineWidth',1.5);
    
    axis tight
    axis square
    
    vline(fixX, 'color', 'k', 'LineStyle',':');
    hline(fixY, 'color', 'k', 'LineStyle',':');
    
    vline(preX, 'color', 'g', 'LineStyle',':');
    hline(preY, 'color', 'g', 'LineStyle',':');
    
    vline(postX, 'color', 'r', 'LineStyle',':');
    hline(postY, 'color', 'r', 'LineStyle',':');
        
    set(gca,'TickDir','out','fontsize',10);
end

% ____________________________________________________________________________ %
%% helper functions
function D = getDist(X,Y)

D = nan(size(X));
D(:,2:end) = sqrt(diff(X,1,2).^2 + diff(Y,1,2).^2);

% ____________________________________________________________________________ %
function mm = minmax(V)
if(~isempty(V))
    mm(1) = nanmin(V(:));
    mm(2) = nanmax(V(:));
else
    mm = [-1 1];
end


