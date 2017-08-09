function  wz_spk_plot_density(spk, timwin, maxResp, EvTm, krnltyp, krnlwidth, maxRast, Rastratio)
% Plot a response estimate based on a spike density function.
%  
% DESCRIPTION 
% This functions plots a response estimate of the firing rate by convolving
% spike events with a spkike density function.
% 
% SYNTAX 
% wz_spk_plot_density(spktimes, tmwin, maxY, EvTm, binwdth, sm, maxRast, Rastratio, clip, trialorder)
%
%   Input:
%         <spktimes>   should be a 2D matrix containing spike timies with rows representing
%                      trials or a structure that has a field 'spiketimes'. 
% 
%         <tmwin>      constrain the shown data to this time window: [min(spk_times(:)) , max(spk_times(:))]
% 
%         <EvTm>       Add data that reflects another event that should be indicated in the data
% 
%         <krnltyp>    Kernel type used as spike density function (default is 'gauss').
% 
%         <krnlwidth>  Kernel width of the spike density function (default is 10 ms).
% 
%         <maxRast>    maximum number of trials represented in the raster plot.
%                      Set it to zero to ommit plotting the raster.
% 
%         <Rastratio>  The portion of the plot that is used for the raster plot.
%                      Set it to zero to ommit plotting the raster.
%
%         <clip>       Discard spike times that occurred after the event specified in <clip>. 
%                      This could be a single number or a vector with one entry per trial.
% 
%         <trialorder> Specify the order of selecting trials. Per default the
%                      trials are ordered according to their row number.
%            
% REFERENCES 
%
% 
%
% ......................................................................... 
% wolf zinke, wolfzinke@gmail.com 
%
% wolf zinke, 14.1.2014

if(exist('timwin','var') == 0 )
    timwin = [];
end

if(exist('krnltyp','var') == 0 )
    krnltyp = 'gauss';
end

if(exist('krnlwidth','var') == 0 )
    krnlwidth = 10;
end

chkpars = 0;
if(~isfield(spk,'spikedensities'))
   chkpars = 1;
elseif(~isempty(krnltyp) && ~isempty(krnlwidth))
    if(~strcmp(krnltyp, spk.kerneltype) || krnlwidth ~= spk.kernelwidth)
        chkpars = 1;
    end
   if(isempty(spk.spikedensities))
       warning('No data available - Nothing to show!');
       return;
   end
end

if(chkpars == 1)
    spk = wz_spk_density(spk, krnltyp, krnlwidth); 
end

if(isempty(timwin) == 1)
    timwin = spk.timewindow;
end
 
nTrials = size(spk.spiketimes,1);

if(exist('maxResp','var') == 0 || isempty(maxResp))
    p = spk.xtime >= timwin(1) &  spk.xtime <= timwin(2);
    maxResp = 1.05 * max(max(spk.meandensity(p) + spk.stedensity(p)));
    if(maxResp == 0)
        maxResp = 1;
    end
end

if(exist('EvTm','var') == 0 || isempty(EvTm))
    EvTm = 0;
end

tidx  = 1:nTrials;

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

do_hold = ishold(gcf);

% plot density with error estimates

m_vals = spk.meandensity;
l_vals = spk.meandensity - spk.stedensity;
h_vals = spk.meandensity + spk.stedensity;
x_vals = spk.xtime;

maxpos = min([find(isnan(m_vals), 1, 'first'), find(isnan(l_vals), 1, 'first'),find(isnan(h_vals), 1, 'first')]);
m_vals(maxpos:end) = [];
l_vals(maxpos:end) = [];
h_vals(maxpos:end) = [];
x_vals(maxpos:end) = []; 

% m_vals = spk.mediandensity;
% % l_vals = spk.mediandensity - spk.IQRerr;
% % h_vals = spk.mediandensity + spk.IQRerr;
% l_vals = spk.lowerdense;
% h_vals = spk.upperdense;

hold on;
if(~isempty(m_vals))
    patch([x_vals, fliplr(x_vals)],[l_vals fliplr(h_vals)],[0.8,0.8,0.8],'EdgeColor',[0.8,0.8,0.8]);
    plot(x_vals, m_vals, 'k', 'linewidth', 2.5);
    plot(x_vals, h_vals, 'Color', [0.5,0.5,0.5]);
    plot(x_vals, l_vals, 'Color', [0.5,0.5,0.5]);
    %vline(EvTm,'color','g', 'linewidth', 2);
end
plot([EvTm EvTm], [0, maxResp],  'b');

%% plot spike raster
if(maxRast > 0 && Rastratio > 0)
    if(~isempty(spk.spiketimes))

        tstep = (Rastratio * maxResp)/(maxRast+2);

        ctrial = maxResp + tstep;
        maxResp = ctrial;
        for(s=1:maxRast)
            if(s > nTrials)
                break;
            end
            ctrial  = ctrial + tstep;
            tpos = spk.trial_order(tidx(s));
            cspikes = spk.spiketimes(tpos,:);
            cspikes(isnan(cspikes)) = [];
            if(~isempty(cspikes))
                plot(cspikes,ctrial,'.k');
            end
            if(length(spk.clip) >= maxRast)
                plot([spk.clip(tpos), spk.clip(tpos)],[ctrial-tstep/2, ctrial+tstep/2],'-r', 'LineWidth',2.5);
            end
        end
        maxResp = maxResp + Rastratio * maxResp;
    end
end

xlim(timwin);
ylim([0, maxResp]);

ylabel('firing rate [spikes/s]');
xlabel('time [ms]');
set(gca,'TickDir','out');

set(gcf, 'Renderer','painters');
set(gcf, 'DoubleBuffer', 'on');

%% plot distribution of clip times
if(isfield(spk,'clip'))
    if(length(spk.clip) == nTrials && nTrials > 0 && spk.clipped == 0)
        [f,xi] = ksdensity(spk.clip,'width', 10); 

        xhist = min(spk.clip)-10 : 20 : max(spk.clip)+10;
        [y,x] = hist(spk.clip,xhist);
        y = (y/max(y)) * 0.1*maxResp;
        f = (f/max(f)) * 0.1*maxResp;
        bar(x,y,1,'r');
        plot(xi,f, 'color', [0.35 0 0] ,'Linewidth',2);

    plot([median(spk.clip),median(spk.clip)], [0, maxResp],  'r');
    %    vline(median(spk.clip),'color','red');
    %    vline(xi(find(f==max(f),1,'first')),'color','r');
    %    vline(mode(spk.clip),'color','r');
    end
end

if(~do_hold)
    hold off;
end



