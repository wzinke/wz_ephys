function SPKWV = wz_spk_get_wave(spikewave, thr, do_plot, ip, ttl, jw)
% Get the widh of a spike waveform.
%
% DESCRIPTION
%   Determine the width of a spike waveform.
%   roughly based on the meanwave function provided by Maik Stüttgen in the MLib toolbox.
%   [http://www.mathworks.com/matlabcentral/fileexchange/37339-mlib-toolbox-for-analyzing-spike-data]
%
% SYNTAX
%
%   [width, Wave] = wz_spk_get_wave(wavemat, do_plot, ip)
%
%   Input:
%
%       wavemat     Vector of mean spike wave or matrix with single spike waves.
%                   Each row corresponds to a spike wave
%
%       thr         specify the threshold in order to select negative (default)
%                   or positive spike amplitude
%
%       do_plot     plot results
%
%       ip          Upsampling factor for the interpolation
%
%       ttl         names used as title in the plot
%
%       jw          plot just the wave, no additional estimates
%
% .........................................................................
% wolf zinke, wolfzinke@gmail.com
%
% $Created : 22-Jan-2015 by wolf zinke
%
% ToDo: - Do Bootstrap
%       - calculate wave forms for specified time periods in order to check stability over time
%       - determine if positive or negative first peak
%       - determine number of phases in the signal
%       - use peak detection to determine trough and preceding and subsequent peaks

% ____________________________________________________________________________ %
%% get information about the current file version
% if(exist('wz_get_git_version','file'))
%     [gitversion, ~, GitInfo] = wz_get_git_version(mfilename('fullpath'));
%     SPKWV.MFileVersion = gitversion;
%     SPKWV.MFileDate    = GitInfo.date;
% end

% ____________________________________________________________________________ %

if(~exist('ip','var') || isempty(ip))
    ip = 100;
end
if(~exist('ip','var'))
    ttl = [];
end

if(nargout == 0)
    do_plot = 1;
end

if(~exist('thr','var') || isempty(thr))
    thr = -1; % assume negative spike waveform as default
end

if(~exist('do_plot','var') || isempty(do_plot))
    do_plot = 0;
end

if(~exist('jw','var') || isempty(jw))
    jw = 0;
end

% ____________________________________________________________________________ %

% subtract median to ensure zero passings
spikewave = double((spikewave - median(spikewave(:))));

if(~isvector(spikewave))
    wavemat   = spikewave;
    spikewave = prctile(wavemat, 50);
    h_val     = prctile(wavemat, 75);
    l_val     = prctile(wavemat, 25);
else
    wavemat = [];
    h_val   = [];
    l_val   = [];
end

med_val = median(wavemat(:));
mad_val = median(abs(wavemat(:) - med_val));

% ____________________________________________________________________________ %
%% interpolate time course
smpltime = 1:length(spikewave);
iptime   = linspace(1,length(spikewave),ip*smpltime(end)) ;

orig_ipwave = spline(smpltime,spikewave,iptime);

if(~isempty(wavemat))
    h_val   = spline(smpltime, h_val, iptime);
    l_val   = spline(smpltime, l_val, iptime);
end
xvec = 1:length(orig_ipwave);

% invert waveform if needed
if(thr > 0)
    ipwave = -1 * orig_ipwave;
    tmphv = h_val;
    h_val  = -1 * l_val;
    l_val  = -1 * tmphv;
    spksgn = -1;
else
    spksgn = 1;
    ipwave = orig_ipwave;
end

% ____________________________________________________________________________ %
%% get spike waveform parameters

% get trough
trough = find(ipwave == min(ipwave) & xvec < xvec(end)/2,1,'first');
if(isempty(trough))
    trough = NaN;
    troughval  = NaN;
else
    troughval  = ipwave(trough);
end
halftrough = troughval/2;

% get max
if(~isnan(trough))
    maxpos = find(ipwave == max(ipwave(trough:end)) & xvec > trough,1,'first');
else
    maxpos = [];
end

if(isempty(maxpos))
    maxpos = NaN;
    maxval = NaN;
else
    maxval = ipwave(maxpos);
end

if(~isnan(trough))
    initmax = find(ipwave == max(ipwave(1:trough)) & xvec < trough,1,'first');
else
    initmax = [];
end

if(isempty(initmax))
    initmax = NaN;
    initmaxvl  = NaN;
else
    initmaxvl  = ipwave(initmax);
end

% get half peak passings
if(~isnan(trough))
    halftrough_down = find(ipwave >= halftrough & xvec < trough, 1,'last');
else
    halftrough_down = [];
end

if(isempty(halftrough_down))
    halftrough_down = NaN;
end

if(~isnan(trough))
    halftrough_up = find(ipwave >= halftrough & xvec > trough, 1,'first');
else
    halftrough_up = [];
end

if(isempty(halftrough_up))
    halftrough_up = NaN;
end

% determine zero passings
if(~isnan(halftrough_down))
    Pre_zeropass  = find(ipwave >= 0 & xvec < halftrough_down, 1,'last');
else
    Pre_zeropass = [];
end

if(isempty(Pre_zeropass))
    Pre_zeropass = NaN;
end

if(~isnan(halftrough_down))
    Post_zeropass = find(ipwave >= 0 & xvec > halftrough_up,   1,'first');
else
    Post_zeropass = [];
end

if(isempty(Post_zeropass))
    Post_zeropass = NaN;
end

% calculate area beluv curve
if(isnan(Pre_zeropass) || isnan(Post_zeropass))
    spkarea = NaN;
else
    spkarea = sum(abs(ipwave(Pre_zeropass:Post_zeropass)));
end
area_abs =  spkarea / ip;
area_rel = (spkarea / sum(abs(ipwave))) / ip;

% calculate width and duration
spike_width    = (halftrough_up - halftrough_down) / ip;
spike_duration = (Post_zeropass - Pre_zeropass)    / ip;

% check if there might have been an inverted spike amplitude
spkinv = mean(ipwave(1:10)) < (-1*troughval);
if(spkinv == 0)
    spkinv = -1;
end

% ____________________________________________________________________________ %
if(nargout == 1)
    SPKWV.wavemat   = wavemat;
    SPKWV.threshold = thr;
    SPKWV.spikewave = orig_ipwave;
    SPKWV.IQR       = spksgn * [h_val, l_val];

    SPKWV.MAD        = mad_val;
    SPKWV.RNG        = spksgn * median([h_val - l_val]);
    SPKWV.peak       = spksgn * troughval;
    SPKWV.minmax     = maxval - troughval;
    SPKWV.peakdist   = maxpos - trough;
    SPKWV.initmax    = initmaxvl - troughval;
    SPKWV.initdist   = trough - initmax;
    SPKWV.SNR        = abs(troughval-med_val) / mad_val;
    SPKWV.SNR2       = abs(troughval-med_val) / SPKWV.RNG;
    SPKWV.troughbin  = trough ./ ip;
    SPKWV.halftrough = [halftrough_down, halftrough_up] ./ ip;
    SPKWV.zeropass   = [Pre_zeropass Post_zeropass]./ip;
    SPKWV.width      = spike_width/ip;
    SPKWV.duration   = spike_duration/ip;
    SPKWV.area       = [area_abs, area_rel]./ip;

    if(thr < 0)
        SPKWV.TreshSign  = 'neg';
    else
        SPKWV.TreshSign  = 'pos';
    end

    if((spkinv * sign(thr)) < 0)
        SPKWV.SpikeSign  = 'neg';
    else
        SPKWV.SpikeSign  = 'pos';
    end
end

% ____________________________________________________________________________ %
if(do_plot)
    hold on;

    h_val = spksgn * h_val;
    l_val = spksgn * l_val;

    if(~isempty(wavemat))
        patch([xvec, fliplr(xvec)],[h_val fliplr(l_val)],[0.8,0.7,0.7], ...
                  'FaceColor',[0.8,0.8,0.8], 'EdgeColor',[0.8,0.8,0.8]);

        Yrng = [min(l_val) max(h_val)];
    else
        Yrng = [min(orig_ipwave) max(orig_ipwave)];
    end

    plot(xvec, orig_ipwave,'k','LineWidth',2);

    hline(0,'LineWidth',0.8);

    if(jw == 0)
        vline(trough,'color','r','LineWidth',2);
        vline([halftrough_down halftrough_up],'color','r');
        vline(maxpos,'color','b','LineWidth',2);
        vline(initmax,'color','b','LineWidth',2);

        plot([halftrough_down, halftrough_up],[halftrough halftrough],'g','LineWidth',2);
    end
    
    if(mod(SPKWV.threshold,1) ~= 0)
        hline(SPKWV.threshold, 'color', 'g','LineStyle',':');
    end
    
    ylim([sort(Yrng) + 0.05*range(Yrng) * [-1 1]]);
    xlim([xvec(1) xvec(end)]);

    title(ttl);
    axis tight
end
