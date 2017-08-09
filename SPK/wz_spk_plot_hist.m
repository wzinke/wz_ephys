function wz_spk_plot_hist(spktimes, tmwin, maxY, EvTm, binwdth, sm, maxRast, Rastratio, clip, trialorder)
% Plot a peri event histogram..
%
% DESCRIPTION
% This functions creates a histogram for the spike occurrences and adds the
% corresponding raster plot.
%
% SYNTAX
% wz_spk_plot_hist(spktimes, tmwin, maxY, EvTm, binwdth, sm, maxRast, Rastratio, clip, trialorder)
%
%   Input:
%         <spktimes>   should be a 2D matrix containing spike timies with rows representing
%                      trials or a structure that has a field 'spiketimes'.
%
%         <tmwin>      constrain the shown data to this time window: [min(spk_times(:)) , max(spk_times(:))]
%
%         <EvTm>       Add data that reflects another event that should be indicated in the data
%
%         <binwdth>    width of the histogram bins.
%
%         <sm>         overlay a smoothed response estimate based on spike convolution with a Gauss function.
%
%         <maxRast>    maximum number of trials represented in the raster plot.
%                      Set it to zero to omit plotting the raster.
%
%         <Rastratio>  The portion of the plot that is used for the raster plot.
%                      Set it to zero to omit plotting the raster.
%
%         <clip>       Discard spike times that occurred after the event specified in <clip>.
%                      This could be a single number or a vector with one entry per trial.
%
%         <trialorder> Specify the order of selecting trials. Per default the
%                      trials are ordered according to their row number.
%
% REFERENCES
%
% .........................................................................
% wolf zinke, wolfzinke@gmail.com
%
% wolf zinke, 14.1.2014
%
%  ToDo: - Make raster plot an extra subplot (optional)
%

%% set default values
if(exist('tmwin','var') == 0 )
    tmwin = [];
end

% make it flexible by also allow structs as input and account for my naming inconsistencies
if(isstruct(spktimes))
    spkstruct = spktimes;
    spktm_alias = {'spktimes', 'spk_times', 'spike_times', 'spiketimes'};
    for(i=1:length(spktm_alias))
        if(isfield(spktimes,spktm_alias(i)))
            spktimes = spktimes.(char(spktm_alias(i)));
        end
    end
else
    spkstruct = [];
end

nTrials = size(spktimes,1);

if(~exist('tmwin','var') || isempty(tmwin))
   tmwin = [floor(min(spktimes(:))) , ceil(max(spktimes(:)))];
else
    spktimes = wz_cropNaN(spktimes, tmwin(1), tmwin(2), NaN);
end

if(~exist('binwdth','var') || isempty(binwdth))
    binwdth = 10;
end

if(exist('EvTm','var') == 0 || isempty(EvTm) == 1)
    EvTm = 0;
end

if(exist('sm','var') == 0 || isempty(sm) == 1)
    sm = 0;
end


tidx = 1:nTrials;

if(exist('maxRast','var') == 0 || isempty(maxRast))
    maxRast = nTrials;
else
    if(nTrials > maxRast)
        tidx  = sort(randperm(nTrials,maxRast));
    end
end

if(exist('Rastratio','var') == 0 || isempty(Rastratio))
    Rastratio = 1/3;
end

%% time vector
xhist = tmwin(1)-binwdth/2 : binwdth : tmwin(2)+binwdth/2;
xvec  = xhist + binwdth/2;


%% clip spike times
if(~exist('clip','var'))
    if(isfield(spkstruct,'clip'))
        clip = spkstruct.clip;
    else
        clip = [];
    end
end

trialcnt = nTrials;
if(~isempty(clip))
    [spktimes, trialorder] = wz_spk_clip(spktimes,clip);

    if(length(clip) == nTrials)
        trialcnt = bsxfun(@le, xvec, clip);
        trialcnt = sum(trialcnt);
    end
end

if(exist('trialorder','var') == 0 || isempty(trialorder))
    trialorder = 1:nTrials;
    if(~isempty(spkstruct))
        if(isfield(spkstruct, 'trial_order'))
            trialorder = spkstruct.trial_order;
        end
    end
end

%% get spike histogram
spkhist = histc(spktimes(:), xhist);
spkhist = (spkhist ./ trialcnt(:)) .* (1000/binwdth); % first term gets the trial average, second one converts spikes/s.

% get scale
if(exist('maxY','var') == 0 || isempty(maxY) == 1)
    maxY = 1.05 * max(spkhist(:));
    if(maxY == 0)
        maxY = 1;
    end
end

%% plot histogram
bar(xvec, spkhist,'k','EdgeColor', 'k', 'BarWidth', 1);
hold on;

if(sm == 1)
%     [f,xi] = ksdensity(spktimes(:),'bandwidth',binwdth);
%     f = f * sum(isnan(spktimes(:)))*nTrials;
    spk = wz_spk_density(spktimes, 'gauss', binwdth);

    sm_ts = spk.meandensity;

    if(length(clip) == nTrials)
        trialcnt = bsxfun(@le, spk.xtime, clip);
        trialcnt = sum(trialcnt);

        sm_ts = (spk.meandensity .* nTrials) ./ trialcnt;
    end

    plot(spk.xtime, sm_ts, 'color', 'b' ,'Linewidth',2);
end

vline(median(EvTm(:)),'color','g', 'linewidth', 2);

%% plot spike raster
if(maxRast > 0 && Rastratio > 0)
    if(~isempty(spktimes))

        tstep = (Rastratio * maxY)/(maxRast+2);

        ctrial = maxY + tstep;
        maxY = ctrial;
        for(s=1:maxRast)
            if(s > nTrials)
                break;
            end
            ctrial  = ctrial + tstep;
            tpos = trialorder(tidx(s));
            cspikes = spktimes(tpos,:);
            cspikes(isnan(cspikes)) = [];
            if(~isempty(cspikes))
                plot(cspikes,ctrial,'.k');
            end
            if(length(clip) == nTrials)
                plot([clip(tpos), clip(tpos)],[ctrial-tstep/2, ctrial+tstep/2],'-r', 'LineWidth',2.5);
            end
        end
        maxY = maxY + Rastratio * maxY;
    end
end

xlim(tmwin);
ylim([0, maxY]);

ylabel('firing rate [spikes/s]');
xlabel('time [ms]');
set(gca,'TickDir','out');

set(gcf, 'Renderer','painters')
set(gcf, 'DoubleBuffer', 'on');

%% plot distribution of clip times
if(length(clip) == nTrials)
    [f,xi] = ksdensity(clip,'width', 10);

    xhist = min(clip)-10 : 10 : max(clip)+10;
    [y,x] = hist(clip,xhist);
    y = (y/max(y)) * 0.1*maxY;
    f = (f/max(f)) * 0.1*maxY;
    bar(x,y,1,'r');
    plot(xi,f, 'color', [1, 0.75,0] ,'Linewidth',2);

    vline(median(clip(:)),'color','red');
%    vline(xi(find(f==max(f),1,'first')),'color','r');
%    vline(mode(spk.clip),'color','r');
end

hold off;

