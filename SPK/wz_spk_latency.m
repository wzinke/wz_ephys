function [lat, latObj] = wz_spk_latency(spks, latwin, spontwin, nBoot, smplsz, krnl, kw, doplot, nfo)
% Latency estimate that defines response onset as outlier compared to the
% base line activity. Detection of outliers utilizes the modified z-score
% as robust measure, which itself is derived from the median absoulute
% deviation (mad).
%
% wolf zinke, July 2014
%
% based on:
%
% B Letham & T Raij (2011)
% Statistically robust measurement of evoked response onset latencies
% Journal of neuroscience methods, 194, 374–379
%
% Leys, C., et al. (2013),
% Detecting outliers: Do not use standard deviation around the mean, use absolute deviation around the median
% Journal of Experimental Social Psychology, 49(4), 764–766.
%
% Boris Iglewicz & David Caster Hoaglin (1993)
% How to Detect and Handle Outliers
% ASQC Quality Press
%
% and half maximum method as described roughly by:
%
% HS Friedman & CE Priebe (1998)
% Estimating stimulus response latency
% Journal of Neuroscience Methods, 83, 185–194
%
% TJ Gawne, TW Kjaer, & BJ Richmond (1996)
% Latency: Another Potential Code for Feature Binding in Striate Cortex
% Journal of Neurophysiology, 76(2), 1356-1360
%
% see also:
% http://www.mathworks.com/matlabcentral/fileexchange/28501-tests-to-identify-outliers-in-data-series
%
% ToDo: This function is not completely checked and has to be considered buggy.

% bootspont =  1;  % define threshold based on bootstrap samples instead of the overall sample
mincons   = 20;  % activity must exceed threshold for at least <mincons> time bins to be considered as response
plotwin   = [-50 250];

%% define default arguments
if(~exist('latwin','var') || isempty(latwin))
    latwin = [30, 200];
end

if(~exist('spontwin','var') || isempty(spontwin))
    spontwin = [-250, 0];
end

if(~exist('nBoot','var') || isempty(nBoot))
    nBoot = 1;
end

if(~exist('smplsz','var') || isempty(smplsz))
    smplsz = size(spks,1);
end

if(~exist('krnl','var') || isempty(krnl))
    krnl = 'gauss';
end

if(~exist('kw','var') || isempty(kw))
    kw = 10;
end

if(~exist('doplot','var') || isempty(doplot))
    doplot = 1;
end

if(~exist('nfo','var'))
    nfo    = [];
elseif(~isempty(nfo))
    doplot = 1;
end

%% get spike density function
spks = wz_cropNaN(spks, min([spontwin(1),latwin(1)]) - 50,  max([spontwin(2),latwin(2)]) + 500);

respest = wz_spk_density(spks, krnl, kw);

% get response estimate
latwinpos = find(respest.xtime >=  latwin(1) & respest.xtime  <  latwin(2));
ppos      =      respest.xtime >= plotwin(1) & respest.xtime <= plotwin(2);
px        =      respest.xtime(ppos);

meanraw  = nanmean( respest.spikedensities);
trialmat = isfinite(respest.spikedensities);
trialcnt = sum(trialmat);
errraw   = nanstd(respest.spikedensities) ./ sqrt(trialcnt);

%% get estimate of spontaneous activity
spontpos  = respest.xtime >= spontwin(1) & respest.xtime < spontwin(2);
% spontrate = 1000 * (nansum(respest.spikedensities(:, spontpos),2) ./ diff(spontwin));
spontrate = nanmean(respest.spikedensities(:, spontpos),2);

%spontrate(isnan(spontrate)) = [];
valid_spont   = isfinite(spontrate);
median_spont  = median(spontrate(valid_spont));
mad_spont     = median(abs(median_spont- spontrate(valid_spont)));  % mad(spontrate(valid_spont),2);
% iqr_spont     = iqr(spontrate(valid_spont));
% P25_spont     = prctile(spontrate(valid_spont),25);
% P75_spont     = prctile(spontrate(valid_spont),75);

sm_spont = prctile(meanraw(spontpos),[1 50 99]);

%% get peak response
maxresp   = max(meanraw(latwinpos));
% maxtm     = max(meanraw(latwinpos));

peakrng = 0.1*(maxresp-median_spont);  % a peak has to be distinct by at least 10% of the range of max response and base line
[ptm, ttm , prt, trt] = peakdet(meanraw, abs(peakrng/4), respest.xtime);

peakpos = find(ptm >= latwin(1) & ptm <= latwin(2));

if(~isempty(peakpos))
    first_peak = ptm(peakpos(1));
    mxpos      = peakpos(prt(peakpos) == max(prt(peakpos)));
    max_peak   = ptm(mxpos(1));
    peakresp   = [prt(peakpos(1)), prt(mxpos(1))];
    
    tp = find(ttm > latwin(1) - 10 & ttm < first_peak);
    if(~isempty(tp))
        first_trough_time = ttm(tp(1));
        first_trough_resp = trt(tp(1));
    else
        first_trough_time = NaN;
        first_trough_resp = NaN;
    end
    
    tp = find(ttm > first_peak & ttm < max_peak(1));
    if(~isempty(tp))
        max_trough_time = ttm(tp(find(trt(tp) == max(trt(tp)) & tp,1,'first')));
        max_trough_resp = trt(tp(find(trt(tp) == max(trt(tp)) & tp,1,'first')));
    else
        max_trough_time = NaN;
        max_trough_resp = NaN;
    end
else
    first_peak = NaN;
    max_peak   = NaN;
    peakresp   = [NaN, NaN];
    
    tp = find(ttm > latwin(1) - 10 & ttm <= latwin(2));
    if(~isempty(tp))
        first_trough_time = ttm(tp(1));
        first_trough_resp = trt(tp(1));
        max_trough_time   = ttm(tp(trt(tp) == max(trt(tp))));
        max_trough_resp   = trt(tp(trt(tp) == max(trt(tp))));
    else
        first_trough_time = NaN;
        first_trough_resp = NaN;
        max_trough_time   = NaN;
        max_trough_resp   = NaN;
    end
end

peakpos = find(ptm >= latwin(1));

if(~isempty(peakpos))
    max_peakresp = ptm(peakpos(prt(peakpos) == max(prt(peakpos))));
else
    max_peakresp = NaN;
end

if(max_peakresp > max_peak(1))
    warning(['Maximal response peak found after latency search window at ',int2str(max_peakresp),' ms!'])
end

%% create a matrix with the bootstrap samples
AllSmp = randi(size(spks,1),   nBoot, smplsz);
if(smplsz == size(spks,1))
    AllSmp(1,:) = 1:smplsz; % force the first value be the raw value
end

%% calculate bootstrapped estimates for the spontaneous activity
boot_median_spont = nanmedian(spontrate(AllSmp),2);
% boot_mad_spont    = mad(spontrate(AllSmp),1,2);

%% get bootstrap estimates for the responses
boot_resp = nan(nBoot, length(respest.xtime));
boot_first_peak        = nan(nBoot, 1);
boot_max_peak          = nan(nBoot, 1);
boot_peakresp          = nan(nBoot, 2);
boot_first_trough_time = nan(nBoot, 1);
boot_max_trough_time   = nan(nBoot, 1);

lw1 = latwin(1);
lw2 = latwin(2);

spkdense = respest.spikedensities;
xtime    = respest.xtime;

for(b=1:nBoot)
    boot_resp(b,:) = nanmean(spkdense(AllSmp(b,:),:));
    
    [ptm, ttm, prt, trt] = peakdet(boot_resp(b,:), abs(peakrng), xtime); % a peak has to be distinct by at least 10% of the range of max response and base line
    
    peakpos = find(ptm >= lw1 & ptm <= lw2);
    
    if(~isempty(peakpos))
        boot_first_peak(b) = ptm(peakpos(1));
        mxpos  = peakpos(find(prt(peakpos) == max(prt(peakpos)),1,'first'));
        boot_max_peak(b)   = ptm(mxpos);
        boot_peakresp(b,:) = [prt(peakpos(1)), prt(mxpos)];
        
        tp = find(ttm > lw1 - 10 & ttm < first_peak);
        if(~isempty(tp))
            boot_first_trough_time(b) = ttm(tp(1));
        end
        
        tp = find(ttm > first_peak & ttm < max_peak(1));
        if(~isempty(tp))
            boot_max_trough_time(b) = ttm(tp(find(trt(tp) == max(trt(tp)) & tp,1,'first')));
        end
    else
        tp = find(ttm > lw1 & ttm < lw2);
        if(~isempty(tp))
            boot_first_trough_time(b) = ttm(tp(1));
            boot_max_trough_time(b)   = ttm(tp(find(trt(tp) == max(trt(tp)) & tp,1,'first')));
        end
    end
end

boot_sm_spont = prctile(boot_resp(:,spontpos), [5 50 95],2);

%% identify super threshold activity
% if(bootspont == 1)
% exceed mad criterion
%    thrvec  = boot_median_spont + 2.3*boot_mad_spont;
thrvec  = boot_sm_spont(:,3);
lat_mat = bsxfun(@gt, boot_resp(:,latwinpos), thrvec(:));

% half max
hfmx     = boot_median_spont + 0.5*(max(boot_resp(:,latwinpos),[],2)-boot_median_spont);
no_resp  = max(boot_resp(:,latwinpos),[],2) > thrvec;
hfmx_mat = bsxfun(@gt, boot_resp(:,latwinpos), hfmx(:));
hfmx_mat = bsxfun(@times, hfmx_mat, no_resp(:));

% half first peak
hfpk     = boot_median_spont + 0.5*(boot_peakresp(:,1)-boot_median_spont);
no_resp  = boot_peakresp(:,1) > thrvec;
hfpk_mat = bsxfun(@gt, boot_resp(:,latwinpos), hfpk(:));
hfpk_mat = bsxfun(@times, hfpk_mat, no_resp(:));

% half max peak
hfmk     = boot_median_spont + 0.5*(boot_peakresp(:,2)-boot_median_spont);
no_resp  = boot_peakresp(:,2) > thrvec;
hfmk_mat = bsxfun(@gt, boot_resp(:,latwinpos), hfmk(:));
hfmk_mat = bsxfun(@times, hfmk_mat, no_resp(:));

% Pre-response supression
thrvec  = boot_sm_spont(:,1);
pep_mat = bsxfun(@lt, boot_resp(:,latwinpos), thrvec(:));

% else
%     % exceed mad criterion
%     lat_mat =  boot_resp(:,latwinpos) > median_spont + 2.3 * mad_spont;
%
%     % half max
%     hfmx_mat = boot_resp(:,latwinpos) > median_spont + 0.5*(maxresp-median_spont);
%
%     % half first peak
%     hfpk_mat = boot_resp(:,latwinpos) > median_spont + 0.5 * (peakresp(1) - median_spont);
%
%     % half max peak
%     hfmk_mat = boot_resp(:,latwinpos) > median_spont + 0.5 * (peakresp(2) - median_spont);
% end

% % multiply with lat_mat to ensure a clear response increase
% hfmx_mat = hfmx_mat.*lat_mat;
% hfpk_mat = hfpk_mat.*lat_mat;
% hfmk_mat = hfmk_mat.*lat_mat;

%% get the first occurences of super threshold activity
boot_lat   = nan(1,nBoot);
boot_hfmx  = nan(1,nBoot);
boot_hfpk  = nan(1,nBoot);
boot_hfmk  = nan(1,nBoot);
boot_pep   = nan(1,nBoot);

boot_rise   = nan(1,nBoot);
boot_decay  = nan(1,nBoot);

xtm = respest.xtime(latwinpos);
valpat = ones(1,mincons);

for(b=1:nBoot)
    p = strfind(lat_mat(b,:), valpat);
    if(~isempty(p))
        boot_lat(b) = xtm(p(1));
        
        if(p(1) > 1)
            pp = find(diff(boot_resp(b, latwinpos(1:p(1))),1)<0,1,'last') + 1;
            if(~isempty(pp))
                boot_rise(b) = xtm(pp);
            end
        end
    end
    
    p = strfind(hfmx_mat(b,:), valpat);
    if(~isempty(p))
        boot_hfmx(b) = xtm(p(1));
    end
    
    p = strfind(hfpk_mat(b,:), valpat);
    if(~isempty(p))
        boot_hfpk(b) = xtm(p(1));
    end
    
    p = strfind(hfmk_mat(b,:), valpat);
    if(~isempty(p))
        boot_hfmk(b) = xtm(p(1));
    end
    
    p = strfind(pep_mat(b,:), ones(1,5)); % use a short period here
    if(~isempty(p))
        if(xtm(p(1)) < boot_lat(b))
            boot_pep(b)   = xtm(p(1));
            if(p(1) > 1)
                pp = find(diff(boot_resp(b, latwinpos(1:p(1))),1)>0,1,'last') + 1;
                if(~isempty(pp))
                    boot_decay(b) = xtm(pp);
                end
            end
        end
    end
end

%% define output
lat = median(boot_hfmx(isfinite(boot_hfmx)));

latObj.nBoot = nBoot;

latObj.mean_resp     = meanraw;
latObj.mean_resp_ste = errraw;

latObj.spontrate = median_spont;
latObj.spontMAD  = mad_spont;

latObj.median_sm_spont = sm_spont(2);
latObj.CI99_sm_spont   = sm_spont([1,3]);

latObj.first_peak   = first_peak;
latObj.max_peak     = max_peak(1);
latObj.peakresp     = peakresp;
latObj.max_peakresp = max_peakresp;

latObj.first_trough_time = first_trough_time;
latObj.first_trough_resp = first_trough_resp;
latObj.max_trough_time   = max_trough_time;
latObj.max_trough_resp   = max_trough_resp;

latObj.halfmax       = median( boot_hfmx(isfinite(boot_hfmx)));
latObj.halfmaxCI95   = prctile(boot_hfmx(isfinite(boot_hfmx)),[ 5 95]);
latObj.halfmaxCI75   = prctile(boot_hfmx(isfinite(boot_hfmx)),[25 75]);
latObj.halfmax_acc   = 100 * sum(isfinite(boot_hfmx)) / nBoot;
latObj.boot_halfmax  = boot_hfmx;

latObj.MADcrit       = median( boot_lat(isfinite(boot_lat)));
latObj.MADcritCI95   = prctile(boot_lat(isfinite(boot_lat)),[ 5 95]);
latObj.MADcritCI75   = prctile(boot_lat(isfinite(boot_lat)),[25 75]);
latObj.MADcrit_acc   = 100 * sum(isfinite(boot_lat)) / nBoot;
latObj.boot_MADcrit  = boot_lat;

latObj.halffrst      = median( boot_hfpk(isfinite(boot_hfpk)));
latObj.halffrstCI95  = prctile(boot_hfpk(isfinite(boot_hfpk)),[ 5 95]);
latObj.halffrstCI75  = prctile(boot_hfpk(isfinite(boot_hfpk)),[25 75]);
latObj.halffrst_acc  = 100 * sum(isfinite(boot_hfpk)) / nBoot;
latObj.boot_halffrst = boot_hfpk;

latObj.Trise         = median( boot_rise(isfinite(boot_rise)));
latObj.TriseCI95     = prctile(boot_rise(isfinite(boot_rise)),[ 5 95]);
latObj.TriseCI75     = prctile(boot_rise(isfinite(boot_rise)),[25 75]);
latObj.Trise_acc     = 100 * sum(isfinite(boot_rise)) / nBoot;
latObj.boot_Trise    = boot_rise;

latObj.halfpeak      = median( boot_hfmk(isfinite(boot_hfmk)));
latObj.halfpeakCI95  = prctile(boot_hfmk(isfinite(boot_hfmk)),[ 5 95]);
latObj.halfpeakCI75  = prctile(boot_hfmk(isfinite(boot_hfmk)),[25 75]);
latObj.halfpeak_acc  = 100 * sum(isfinite(boot_hfmk)) / nBoot;
latObj.boot_halfpeak = boot_hfmk;

latObj.maxpeak       = median( boot_max_peak(isfinite(boot_max_peak)));
latObj.maxpeakCI95   = prctile(boot_max_peak(isfinite(boot_max_peak)),[ 5 95]);
latObj.maxpeakCI75   = prctile(boot_max_peak(isfinite(boot_max_peak)),[25 75]);
latObj.maxpeak_acc   = 100 * sum(isfinite(boot_max_peak)) / nBoot;
latObj.boot_maxpeak  = boot_max_peak;

latObj.frstpeak      = median( boot_first_peak(isfinite(boot_first_peak)));
latObj.frstpeakCI95  = prctile(boot_first_peak(isfinite(boot_first_peak)),[ 5 95]);
latObj.frstpeakCI75  = prctile(boot_first_peak(isfinite(boot_first_peak)),[25 75]);
latObj.frstpeak_acc  = 100 * sum(isfinite(boot_first_peak)) / nBoot;
latObj.boot_frstpeak = boot_first_peak;

latObj.firsttrough      = median( boot_first_trough_time(isfinite(boot_first_trough_time)));
latObj.firsttroughCI95  = prctile(boot_first_trough_time(isfinite(boot_first_trough_time)),[ 5 95]);
latObj.firsttroughCI75  = prctile(boot_first_trough_time(isfinite(boot_first_trough_time)),[25 75]);
latObj.firsttrough_acc  = 100 * sum(isfinite(boot_first_trough_time)) / nBoot;
latObj.boot_firsttrough = boot_first_trough_time;

latObj.maxtrough        = median( boot_max_trough_time(isfinite(boot_max_trough_time)));
latObj.maxtroughCI95    = prctile(boot_max_trough_time(isfinite(boot_max_trough_time)),[ 5 95]);
latObj.maxtroughCI75    = prctile(boot_max_trough_time(isfinite(boot_max_trough_time)),[25 75]);
latObj.maxtrough_acc    = 100 * sum(isfinite(boot_max_trough_time)) / nBoot;
latObj.boot_maxtrough   = boot_max_trough_time;

latObj.PEP      = median( boot_pep(isfinite(boot_pep)));
latObj.PEPCI95  = prctile(boot_pep(isfinite(boot_pep)),[ 5 95]);
latObj.PEPCI75  = prctile(boot_pep(isfinite(boot_pep)),[25 75]);
latObj.PEP_acc  = 100 * sum(isfinite(boot_pep)) / nBoot;
latObj.boot_PEP = boot_pep;

latObj.Tdecay      = median( boot_decay(isfinite(boot_decay)));
latObj.TdecayCI95  = prctile(boot_decay(isfinite(boot_decay)),[ 5 95]);
latObj.TdecayCI75  = prctile(boot_decay(isfinite(boot_decay)),[25 75]);
latObj.Tdecay_acc  = 100 * sum(isfinite(boot_decay)) / nBoot;
latObj.boot_Tdecay = boot_decay;


%% plot
if(doplot == 1)
    figure('Position', [0 0 1200 1000], 'Renderer', 'Painters');
    
    %% raster overview
    subplot('Position',[0.05 0.525  0.40 0.425],'fontsize',12);
    hold on;
    title(nfo);
    xlabel('time [ms]')
    ylabel('trial number')
    set(gca,'TickDir','out');
    
    maxY = respest.nTrials+1;
    
    if(~isempty(spks))
        for(t=1:respest.nTrials)
            cspk = spks(t,:);
            cspk(cspk < plotwin(1) | cspk > plotwin(2) | isnan(cspk)) = [];
            if(~isempty(cspk))
                plot(cspk, t, '.k')
            end
        end
    end
    
    vline(latObj.maxtrough, 'color', 'b', 'LineWidth', 2)
    vline(latObj.MADcrit,   'color', 'm', 'LineWidth', 2)
    vline(latObj.halfmax,   'color', 'r', 'LineWidth', 2)
    vline(latObj.frstpeak,  'color', 'c', 'LineWidth', 2)
    vline(latObj.maxpeak,   'color', 'g', 'LineWidth', 2)
    vline(0, 'color', 'k', 'LineWidth', 2)
    
    ylim([0; maxY]);
    xlim(plotwin);
    
    %% surfed response
    subplot('Position',[0.05 0.05  0.40 0.425],'fontsize',12);
    hold on;
    xlabel('time bin [ms]')
    ylabel('trial number')
    set(gca,'TickDir','out');
    
    % ipk = ones(5);
    ipk = ones(15);
    for(b=1:round(size(ipk,1)))
        ipk(b:end-b,b:end-b) = b;
    end
    
    % ipresp  =conv2(respest.spikedensities,ipk,'same');
    ipresp = conv2(respest.spiketrain,ipk,'same');
    
    pcolor(ipresp);
    shading interp;
    
    vline(find(respest.xtime==0), 'color', 'k', 'LineWidth', 2)
    vline(find(respest.xtime==latObj.MADcrit),   'color', 'm', 'LineWidth', 2)
    vline(find(respest.xtime==latObj.halfmax),   'color', 'r', 'LineWidth', 2)
    vline(find(respest.xtime==latObj.maxpeak),   'color', 'g', 'LineWidth', 2)
    vline(find(respest.xtime==latObj.maxtrough), 'color', 'b', 'LineWidth', 2)
    vline(find(respest.xtime==latObj.frstpeak),  'color', 'c', 'LineWidth', 2)
    
    xlim(find(respest.xtime == plotwin(1) | respest.xtime == plotwin(2)));
    ylim([1; maxY-1]);
    
    %% response estimates without clipping
    subplot('Position',[0.5 0.68 0.48 0.275],'fontsize',12);
    hold on;
    ylabel('activity [spikes/s]')
    xlabel('time [s]')
    set(gca,'TickDir','out');
    
    %     uthr = median_spont + 2.3 * mad_spont;
    %     lthr = median_spont - 2.3 * mad_spont;
    
    patch(plotwin([1 2 2 1])+[-10 10 10 -10], sm_spont([1 1 3 3]),[255 190  56]./255,'EdgeColor',[255 125   0]./255);
    hline(sm_spont(2), 'color', [255 125   0]./255, 'LineWidth', 2)
    
    hline(sm_spont(2) + 0.5*(maxresp-sm_spont(2)), 'color', 'r', 'LineWidth', 2)
    
    patch([px, fliplr(px)],[meanraw(ppos)-errraw(ppos), fliplr(meanraw(ppos)+errraw(ppos))], [109 121 235]./255, 'EdgeColor', [  1   5  45]./255);
    plot(px, meanraw(ppos), 'Color',  [  1   5  45]./255, 'linewidth', 2.5);
    
    vline(latObj.maxtrough, 'color', 'b', 'LineWidth', 2)
    vline(latObj.MADcrit,   'color', 'm', 'LineWidth', 2)
    vline(latObj.halfmax,   'color', 'r', 'LineWidth', 2)
    vline(latObj.frstpeak,  'color', 'c', 'LineWidth', 2)
    vline(latObj.maxpeak,   'color', 'g', 'LineWidth', 2)
    vline(0, 'color', 'k', 'LineWidth', 2)
    
    xlim(plotwin);
    
    ylim([0; max(meanraw+errraw)]);
    
    %% boot median outlier criterion
    subplot('Position',[0.50 0.45 0.14 0.16],'fontsize',12);
    hold on;
    ylabel('count');
    title('mod z threshold');
    set(gca,'TickDir','out');
    
    [t, xh] = hist(latObj.boot_MADcrit,min(latObj.boot_MADcrit):max(latObj.boot_MADcrit));
    
    if(~isempty(xh))
        kbw = min(diff(xh));
        bar(xh,t,1,'k','EdgeColor','k')
        [k, x] = ksdensity(latObj.boot_MADcrit, 'width', kbw);
        plot(x,k*kbw*sum(isfinite(latObj.boot_MADcrit)),'-m','LineWidth',2.5);
        
        vline(latObj.MADcrit,     'color', 'm', 'LineWidth', 2)
        vline(latObj.MADcritCI75, 'color', 'm', 'linestyle', '--','LineWidth', 2)
        vline(latObj.MADcritCI95, 'color', 'm', 'linestyle',  ':','LineWidth', 2)
        
        xlim(prctile(latObj.boot_MADcrit,[5,95])+[-2.5 2.5]); %[xh(1)-0.5, xh(end)+0.5]);
        ylim([0, max([max(t), max(k)])]);
    else
        ylim([0, 1]);
        xlim([0, 1]);
    end
    
    %% boot half max
    subplot('Position',[0.675 0.45 0.14 0.16],'fontsize',12);
    hold on;
    title('half max');
    set(gca,'TickDir','out');
    
    [t, xh] = hist(latObj.boot_halfmax,min(latObj.boot_halfmax):max(latObj.boot_halfmax));
    
    if(~isempty(xh))
        kbw = min(diff(xh));
        bar(xh,t,1,'k','EdgeColor','k')
        [k, x] = ksdensity(latObj.boot_halfmax, 'width', kbw);
        plot(x,k*kbw*sum(isfinite(latObj.boot_halfmax)),'-r','LineWidth',2.5);
        
        vline(latObj.halfmax,     'color', 'r', 'LineWidth', 2)
        vline(latObj.halfmaxCI75, 'color', 'r', 'linestyle', '--','LineWidth', 2)
        vline(latObj.halfmaxCI95, 'color', 'r', 'linestyle',  ':','LineWidth', 2)
        
        xlim(prctile(latObj.boot_halfmax,[5,95])+[-2.5 2.5]); %[xh(1)-0.5, xh(end)+0.5]);
        ylim([0, max([max(t), max(k)])]);
    else
        ylim([0, 1]);
        xlim([0, 1]);
    end
    
    %% boot T-rise
    subplot('Position',[0.85 0.45 0.14 0.16],'fontsize',12);
    hold on;
    title('T-rise');
    set(gca,'TickDir','out');
    
    [t, xh] = hist(latObj.boot_Trise,min(latObj.boot_Trise):max(latObj.boot_Trise));
    
    if(~isempty(xh))
        kbw = min(diff(xh));
        bar(xh,t,1,'k','EdgeColor','k')
        [k, x] = ksdensity(latObj.boot_Trise, 'width', kbw);
        plot(x,k*kbw*sum(isfinite(latObj.boot_Trise)),'-r','LineWidth',2.5);
        
        vline(latObj.Trise,     'color', 'r', 'LineWidth', 2)
        vline(latObj.TriseCI75, 'color', 'r', 'linestyle', '--','LineWidth', 2)
        vline(latObj.TriseCI95, 'color', 'r', 'linestyle',  ':','LineWidth', 2)
        
        xlim(prctile(latObj.boot_Trise,[5,95])+[-2.5 2.5]); %[xh(1)-0.5, xh(end)+0.5]);
        ylim([0, max([max(t), max(k)])]);
    else
        ylim([0, 1]);
        xlim([0, 1]);
    end
    
    %% first peak
    subplot('Position',[0.50 0.25 0.14 0.16],'fontsize',12);
    hold on;
    %    xlabel('time [ms]');
    ylabel('count');
    title('first peak response');
    set(gca,'TickDir','out');
    
    [t, xh] = hist(latObj.boot_frstpeak,min(latObj.boot_frstpeak):max(latObj.boot_frstpeak));
    
    if(~isempty(xh))
        kbw = min(diff(xh));
        bar(xh,t,1,'k','EdgeColor','k')
        [k, x] = ksdensity(latObj.boot_frstpeak, 'width', kbw);
        plot(x,k*kbw*sum(isfinite(latObj.boot_frstpeak)),'-c','LineWidth',2.5);
        
        vline(latObj.frstpeak,     'color', 'c', 'LineWidth', 2)
        vline(latObj.frstpeakCI75, 'color', 'c', 'linestyle', '--','LineWidth', 2)
        vline(latObj.frstpeakCI95, 'color', 'c', 'linestyle',  ':','LineWidth', 2)
        
        xlim(prctile(latObj.boot_frstpeak,[5,95])+[-2.5 2.5]); %[xh(1)-0.5, xh(end)+0.5]);
        ylim([0, max([max(t), max(k)])]);
    else
        ylim([0, 1]);
        xlim([0, 1]);
    end
    
    %% max peak
    subplot('Position',[0.675 0.25 0.14 0.16],'fontsize',12);
    hold on;
    %     ylabel('count');
    %     xlabel('time [ms]');
    title('max peak response');
    set(gca,'TickDir','out');
    
    [t, xh] = hist(latObj.boot_maxpeak,min(latObj.boot_maxpeak):max(latObj.boot_maxpeak));
    
    if(~isempty(xh))
        kbw = min(diff(xh));
        bar(xh,t,1,'k','EdgeColor','k')
        [k, x] = ksdensity(latObj.boot_maxpeak, 'width', kbw);
        plot(x,k*kbw*sum(isfinite(latObj.boot_maxpeak)),'-g','LineWidth',2.5);
        
        vline(latObj.maxpeak,     'color', 'g', 'LineWidth', 2)
        vline(latObj.maxpeakCI75, 'color', 'g', 'linestyle', '--','LineWidth', 2)
        vline(latObj.maxpeakCI95, 'color', 'g', 'linestyle',  ':','LineWidth', 2)
        
        xlim(prctile(latObj.boot_maxpeak,[5,95])+[-2.5 2.5]); %[xh(1)-0.5, xh(end)+0.5]);
        ylim([0, max([max(t), max(k)])]);
    else
        ylim([0, 1]);
        xlim([0, 1]);
    end
    
    %% trough
    subplot('Position',[0.85 0.25 0.14 0.16],'fontsize',12);
    hold on;
    xlabel('time [ms]');
    title('trough response suppression');
    set(gca,'TickDir','out');
    
    [t, xh] = hist(latObj.boot_maxtrough,min(latObj.boot_maxtrough):max(latObj.boot_maxtrough));
    
    if(~isempty(xh))
        kbw = min(diff(xh));
        bar(xh,t,1,'k','EdgeColor','k')
        [k, x] = ksdensity(latObj.boot_maxtrough, 'width', kbw);
        plot(x,k*kbw*sum(isfinite(latObj.boot_maxtrough)),'-b','LineWidth',2.5);
        
        vline(latObj.maxtrough,     'color', 'b', 'LineWidth', 2)
        vline(latObj.maxtroughCI75, 'color', 'b', 'linestyle', '--','LineWidth', 2)
        vline(latObj.maxtroughCI95, 'color', 'b', 'linestyle',  ':','LineWidth', 2)
        
        xlim(prctile(latObj.boot_maxtrough,[5,95])+[-2.5 2.5]); %[xh(1)-0.5, xh(end)+0.5]);
        ylim([0, max([max(t), max(k)])]);
    else
        ylim([0, 1]);
        xlim([0, 1]);
    end
    
    %% PEP
    subplot('Position',[0.50 0.050 0.14 0.16],'fontsize',12);
    hold on;
    ylabel('count');
    xlabel('time [ms]');
    title('PEP');
    set(gca,'TickDir','out');
    
    [t, xh] = hist(latObj.boot_PEP,min(latObj.boot_PEP):max(latObj.boot_PEP));
    
    if(~isempty(xh))
        kbw = min(diff(xh));
        bar(xh,t,1,'k','EdgeColor','k')
        [k, x] = ksdensity(latObj.boot_PEP, 'width', kbw);
        plot(x,k*kbw*sum(isfinite(latObj.boot_PEP)),'color', [0.5 0.5 0.5],'LineWidth',2.5);
        
        vline(latObj.PEP,     'color', [0.5 0.5 0.5], 'LineWidth', 2)
        vline(latObj.PEPCI75, 'color', [0.5 0.5 0.5], 'linestyle', '--','LineWidth', 2)
        vline(latObj.PEPCI95, 'color', [0.5 0.5 0.5], 'linestyle',  ':','LineWidth', 2)
        
        xlim(prctile(latObj.boot_PEP,[5,95])+[-2.5 2.5]); %[xh(1)-0.5, xh(end)+0.5]);
        ylim([0, max([max(t), max(k)])]);
    else
        ylim([0, 1]);
        xlim([0, 1]);
    end
    
    
    %% T-Decay
    subplot('Position',[0.675 0.050 0.14 0.16],'fontsize',12);
    
    hold on;
    xlabel('time [ms]');
    title('T-Decay');
    set(gca,'TickDir','out');
    
    [t, xh] = hist(latObj.boot_Tdecay,min(latObj.boot_Tdecay):max(latObj.boot_Tdecay));
    if(~isempty(xh))
        kbw = min(diff(xh));
        bar(xh,t,1,'k','EdgeColor','k')
        [k, x] = ksdensity(latObj.boot_Tdecay, 'width', kbw);
        plot(x,k*kbw*sum(isfinite(latObj.boot_Tdecay)),'color', [0.5 0.5 0.5],'LineWidth',2.5);
        
        vline(latObj.Tdecay,     'color', [0.5 0.5 0.5], 'LineWidth', 2)
        vline(latObj.TdecayCI75, 'color', [0.5 0.5 0.5], 'linestyle', '--','LineWidth', 2)
        vline(latObj.TdecayCI95, 'color', [0.5 0.5 0.5], 'linestyle',  ':','LineWidth', 2)
        
        xlim(prctile(latObj.boot_Tdecay,[5,95])+[-2.5 2.5]); %[xh(1)-0.5, xh(end)+0.5]);
        ylim([0, max([max(t), max(k)])]);
    else
        ylim([0, 1]);
        xlim([0, 1]);
    end
end




