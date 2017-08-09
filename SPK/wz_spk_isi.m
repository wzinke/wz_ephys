function ISIobj = wz_spk_isi(spktimes, timwin, doplot, nfo)
% Calculation of interspike intervals based on times of spike occurrence.
% Input is expected to be a 2-dimensional matrix where the rows correspond to
% trials.
%
% Input:
%       <spktimes>  - vector or 2D-Matrix containing spike times
%       <timewin>   - specifies a time window for calculating ISIs (default over all spikes);
%                     The window could be defined as single window using a start and end time, or
%                     as a set of windows in a vector or 2D matrix. Given times are including values.
% Output:
%       <ISIobj>    - structure summarizing several 
%
% If a 2D matrix is used missing values should be filled up with NaN's.
% It is assumed, that <spktimes> are used to specify the times of the intervals
% (either preceeding or following the given time).
%
% References:
%
% [-] Moore GP, Perkel DH, and Segundo JP (1966)
%     Statistical Analysis and Functional Interpretation of Neuronal Spike Data.
%     Annu Rev Physiol, 28, 493-522
%
% [-] Perkel DH, Gerstein GL, and Moore GP (1967)
%     Neuronal Spike Trains I: The Single Spike Train.
%     Biophysical J, 7, 391-418
%
% [-] Softky WR & Koch C (1993)
%     The High Irregular Firing of Cortical Cells is Inconsistent with Temporal Integration of Random EPSPs.
%     J Neurosci 13(1), 334-350
%
% [-] Holt GR, Softky WR, Koch C, and Douglas RJ (1996).
%     Comparison of Discharge Variability In Vitro and In Vivo in Cat Visual Neurons.
%     J Neurophysiol, 75(5), 1806-1814
%
% [-] Davies RM, Gerstein GL, and Baker SN (2006).
%     Measurement of Time-Dependent Changes in the Irregularity of Neural Spiking.
%     J Neurophysiol, 96(2), 906-18.
%
% [-] Shinomoto S, Shima K, and Tanji J (2003).
%     Differences in Spiking Patterns among Cortical Neurons.
%     Neural Computation, 15, 2823-2842
%
% [-] Shinomoto S, Miura K, Koyama S (2005).
%     A Measure of Local Variation of Inter-Spike Intervals.
%     Biosystems,79(1-3), 67-72
%
% [-] Nowak LG, Azouz R, Sanchez-Vives MV, Gray CM, & McCormick DA (2003).
%     Electrophysiological classes of cat primary visual cortical neurons 
%     in vivo as revealed by quantitative analyses.
%     J Neurophysiol, 89(3), 1541-1566
%
% wolf zinke, 3.9.2006 (reworked May 2014)
%
% ToDo: - add some bootstrapping approach (across trials?) for FF/CV/LV/IR ...
%       - make option for switching between simple ISI analysis and complete analysis
%       - make a reasonable summary plot

%% preparation
spktimes = squeeze(spktimes);   % reduce dimensions

if (exist('timwin','var') == 0)
     timwin = [];
end

if(~exist('doplot','var') || isempty(doplot))
    doplot = 0;
end

if(nargout == 0)
    doplot = 1;
end

if(~exist('nfo','var'))
    nfo    = [];
elseif(~isempty(nfo))
    doplot = 1;
end

%% calculate ISI for each trial
% extract times from specified windows
if(~isempty(timwin))
    spktimes = wz_cropNaN(spktimes, timwin(1), timwin(2));
else
    timwin = [min(spktimes(:)), max(spktimes(:))];
end

spkcnt = sum(isfinite(spktimes),2);

%% get ISIs - it is just the difference of spike times
[isi_vec, isi] = SPK_get_ISI(spktimes);

%% get an ISI histogram
isi_x    = 0:max(isi(:));
isi_hist = hist(isi_vec, isi_x)./length(isi_vec);

%% assign values to output struct
ISIobj.spkcnt   = spkcnt;

ISIobj.NumZero = sum(isi_vec == 0);

ISIobj.spkrates = 1000 * spkcnt ./ (diff(timwin));  
ISIobj.spkrates_mean  = mean(ISIobj.spkrates);
ISIobj.spkrates_std   = std(ISIobj.spkrates);
ISIobj.spkrates_ste   = ISIobj.spkrates_std / sqrt(length(ISIobj.spkrates));

ISIobj.spkrates_median = median(ISIobj.spkrates);
ISIobj.spkrates_mad    = mad(ISIobj.spkrates, 1);
ISIobj.spkrates_iqr    = iqr(ISIobj.spkrates);
ISIobj.spkrates_P25    = prctile(ISIobj.spkrates, [25, 75]);
ISIobj.spkrates_P95    = prctile(ISIobj.spkrates, [ 5, 95]);

ISIobj.isi      = isi;
ISIobj.isi_vec  = isi_vec;
ISIobj.isi_x    = isi_x;
ISIobj.isi_hist = isi_hist;

ISIobj.trial_isi = nanmean(isi,2);
% ISIobj.isi_trial_trend = smooth(ISIobj.trial_isi, 50/ntrial, 'loess');

%% calculate some ISI based statistics
ISIobj.ISImean = mean(isi_vec);

ISIobj.FF  = SPK_get_Fano(spkcnt);

ISIobj.CV  = SPK_get_ISIcv( isi, 0);
ISIobj.CV2 = SPK_get_ISIcv2(isi, 0);   % Holt et al., Softky & Koch
ISIobj.LV  = SPK_get_ISIlv( isi, 0);   % Shinomoto et al.
ISIobj.IR  = SPK_get_ISIir( isi, 0);   % Davies et al.

ISIobj.trial_CV  = SPK_get_ISIcv( isi, 1);
ISIobj.trial_CV2 = SPK_get_ISIcv2(isi, 1);
ISIobj.trial_LV  = SPK_get_ISIlv( isi, 1);
ISIobj.trial_IR  = SPK_get_ISIir( isi, 1);

if(~isempty(isi_hist))
    ISIobj.RefRat1    = sum(isi_hist(1:2)) / sum(isi_hist);
    ISIobj.RefRat2    = sum(isi_hist(1:2)) / max(isi_hist);
    ISIobj.HalfMaxBin = isi_x(find(isi_hist > max(isi_hist)/2, 1, 'first'));
else
    ISIobj.RefRat1    = NaN;
    ISIobj.RefRat2    = NaN;
    ISIobj.HalfMaxBin = NaN;
end

%% log ISI analysis (according to Nowak et al.)
valISI = isi_vec > 0;

ISIobj.logISI = log(isi_vec);

ISIobj.logISI_mean   = mean(ISIobj.logISI(valISI));
ISIobj.logISI_median = median(ISIobj.logISI(valISI));
%ISIobj.logISI_geom   = geomean(ISIobj.logISI(valISI));
ISIobj.logISI_CV     = std(ISIobj.logISI(valISI)) / mean(ISIobj.logISI(valISI));
ISIobj.logISI_skew   = skewness(ISIobj.logISI(valISI));
ISIobj.logISI_kurtos = kurtosis(ISIobj.logISI(valISI));
ISIobj.logISI_iqr    = iqr(ISIobj.logISI(valISI));

try
[hdt_dip, hdt_p] = HartigansDipSignifTest(ISIobj.logISI(valISI), 1);
    ISIobj.hdt_dip   = hdt_dip;
    ISIobj.hdt_p     = hdt_p;
catch
    ISIobj.hdt_dip   = NaN;
    ISIobj.hdt_p     = NaN;    
end
ISIobj.logISI_BC = wz_bimod_test(ISIobj.logISI(valISI));


%% fit distributions (these are crappy fits, need improvement)
% try
%     gamma_fit = fitdist(isi_vec','Gamma');
%     ISIobj.gamma_fit = gamma_fit;
%     ISIobj.gamma_a   = gamma_fit.a;
%     ISIobj.gamma_b   = gamma_fit.b;
% catch
%     ISIobj.gamma_fit = [];
%     ISIobj.gamma_a   = NaN;
%     ISIobj.gamma_b   = NaN;
% end
% 
% try
%     tmpisi = isi_vec;
%     tmpisi(tmpisi == 0) = 1; % correct, otherwise lognorm fails
%     lognorm_fit  = fitdist(tmpisi','Lognormal');
%     ISIobj.lognorm_fit   = lognorm_fit;
%     ISIobj.lognorm_mu    = lognorm_fit.mu;
%     ISIobj.lognorm_sigma = lognorm_fit.sigma;
% catch
%     ISIobj.lognorm_fit   = [];
%     ISIobj.lognorm_mu    = NaN;
%     ISIobj.lognorm_sigma = NaN;
% end

% %% fit gamma distribution
% funcall = 'wz_gamma_mod(x, amp, kapa, tau, lag)';
% 
% lowbnd = [2*min(isi_hist),-10,0.5,0.5*min(isi_hist),0.2];
% uppbnd = [2*max(isi_hist),10,25,0.5*max(isi_hist),10 ];
% initguess = [0.8*(max(isi_hist)) 1 3 mean(isi_hist) 1];
% 
% 
% sf = fitoptions('Method','NonlinearLeastSquares',...
%     'Startpoint',initguess,...
%     'Lower',lowbnd,...
%     'Upper',uppbnd,...
%     'MaxFunEvals',100000,...
%     'MaxIter',100000, ...
%     'Robust','on');
% 
% ft = fittype(funcall,'options',sf);
% [cfun, gof, output] = fit(isi_x,isi_hist,ft);
% 


if(doplot == 1)
    %% Needs some work-over
    hold on
    title(nfo);
    
    [t, xh] = hist(ISIobj.isi_vec,0:max(ISIobj.isi_vec));
    kbw = min(diff(xh));
    
    bar(xh,t,1,'k','EdgeColor','k')
    % bar(isi_x,isi_hist,1,'k','EdgeColor','w');
    
%     gf = pdf(gamma_fit,isi_x);
%     plot(isi_x,gf*length(ISIobj.isi_vec),'LineWidth',2,'color','blue')
%     
%     lf = pdf(lognorm_fit,isi_x);
%     plot(isi_x,lf*length(ISIobj.isi_vec),'LineWidth',2,'color','red')
    
    [k, x] = ksdensity(ISIobj.isi_vec, 0:0.1:max(ISIobj.isi_vec), 'width', kbw);
    plot(x,k*kbw*length(ISIobj.isi_vec),'color', 'g','LineWidth',2.5);
    
    vline(mean(isi_vec),  'color','green');
    vline(median(isi_vec),'color','cyan');

    xlim([0,prctile(isi_vec,85)]);
end

