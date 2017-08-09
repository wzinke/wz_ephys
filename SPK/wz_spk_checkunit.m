function Ustats = wz_spk_checkunit(spk, spontwin, nfo)
% use spontaneous activity to asses recording stability
% and whether the spiking activity likely is from a single
% unit or from multi unit.
%
% also check;
%     DN Hill, SB Mehta, and D Kleinfeld
%     Quality Metrics to Accompany Spike Sorting of Extracellular Signals
%     Journal of Neuroscience, 31(24): 8699-8705
%     doi: 10.1523/JNEUROSCI.0971-11.2011

% nBoot  = 25000;  % we have largge sample sizes, thus need large bootstrap samples
% BootSZ = 200;    % in case to run m by n bootstrap
olblk  = 20;

%binwd  = 50;
peakdelta = 0.1;
smtrials  = 75; % number of trials used for the Loess smoothing
minbblock = 50; % minimum block to be accounted as recording instability

%% check data input
if(~exist('spontwin','var') || isempty(spontwin))
    spontwin = [-500; 0];
end

if(~exist('nfo','var'))
    nfo = [];
end

if(~isstruct(spk) )
    if(length(unique(spk(~isnan(spk)))) > 2)
        spk = wz_get_SPKobj(spk,spontwin);
    end
end

if(spk.timewindow(1) > spontwin(1))
    warning('Selected Time Window does not match time window of spike times!');
end

%% get spiketimes only for the spontaneous activity
sptspk = spk.spiketimes;

% sptspk(sptspk>spontwin(2) | sptspk < spontwin(1)) = NaN;
sptspk = wz_cropNaN(sptspk, spontwin(1), spontwin(2));
spkcnt = sum(isfinite(sptspk),2);

%spbins = spontwin(1) : binwd : spontwin(2);

spontspk  = nan(spk.nTrials,max(spkcnt));
% spontbins = nan(spk.nTrials,length(spbins));
for(t=1:spk.nTrials)
    p = isfinite(sptspk(t,:));
    spontspk(t,1:sum(p)) = sptspk(t,p);

%     spontbins(t,:) = histc(sptspk(t,p),spbins);
end

sprt = 1000*(spkcnt/diff(spontwin));  % spontaneous firing rate for each trial
sm_sprt = smooth(sprt, smtrials/spk.nTrials,'loess');

% spt_mad = mad(sprt,1);
% spt_med = median(sprt);
spt_rng = [min(sprt) , max(sprt)];

if(diff(spt_rng) == 0)
    spt_rng = [0 1];
end

Ustats.sprt     = sprt;
Ustats.sm_sprt  = sm_sprt;
Ustats.spontspk = spontspk;
Ustats.spontwin = spontwin;

% binrate = 1000*(mean(spontbins,2)/binwd);
% sprt = binrate;

% %% bootstrap the spontaneous rate
% ctrials      = randi(spk.nTrials,spk.nTrials,nBoot); %wz_sample(trialvec, num_trial, 1, rep);
% %ctrials      = randi(BootSZ,spk.nTrials,nBoot); %wz_sample(trialvec, num_trial, 1, rep);
% ctrials(:,1) = 1:spk.nTrials;
%
% bootsmpls = sprt(ctrials);
% bootmean  = mean(bootsmpls);

%% get interspike interval information
ISIobj = wz_spk_isi(spontspk, spontwin);

Ustats.isi             = ISIobj.isi;
Ustats.isi_vec         = ISIobj.isi_vec;
Ustats.isi_x           = ISIobj.isi_x;
Ustats.isi_hist        = ISIobj.isi_hist;
% Ustats.isi_trial_trend = ISIobj.isi_trial_trend;
Ustats.isi_IR          = ISIobj.IR;
Ustats.isi_CV          = ISIobj.CV;
Ustats.isi_CV2         = ISIobj.CV2;
Ustats.isi_FF          = ISIobj.FF;

Ustats.isi1stBins2Peak = ISIobj.RefRat2;
Ustats.isi1stBins2all  = ISIobj.RefRat1;
Ustats.ISIHalfMaxBin   = ISIobj.HalfMaxBin;
Ustats.trial_CV        = ISIobj.trial_CV;
Ustats.trial_CV2       = ISIobj.trial_CV2;
Ustats.trial_LV        = ISIobj.trial_LV;
Ustats.trial_IR        = ISIobj.trial_IR;

ISIobj.trial_isi(isnan(ISIobj.trial_isi)) = -1;
sm_isi = smooth(ISIobj.trial_isi, smtrials/spk.nTrials,'loess');


%% plot raster
figure('Position', [0 0 1600 1000], 'Renderer', 'Painters');

subplot('Position',[0.05 0.075  0.25 0.9]);
hold on;

if(~isempty(spontspk))
    for(t=1:spk.nTrials)
        plot(spontspk(t,:),t,'.k')
    end
end
ylim([0; spk.nTrials+1]);
xlim(spontwin);
set(gca,'TickDir','out');
ylabel('trial number');
xlabel('time [ms]')

%% plot histograms
% distribution of raw spontaneous rates
subplot('Position',[0.545 0.71  0.225 0.265]);
hold on;
[t, xh] = hist(sprt,25);

bar(xh,t,1,'k','EdgeColor','k');
[k, x] = ksdensity(sprt, 'width', max(diff(xh)));
plot(x,k*min(diff(xh))*length(sprt),'-r','LineWidth',2.5);

% vline(spt_med,'color','red','LineWidth',2.5);
% vline(spt_med + spt_mad,'color','red','LineStyle','--','LineWidth',2.5);
% vline(spt_med - spt_mad,'color','red','LineStyle','--','LineWidth',2.5);
% vline(uthr,'color','red','LineStyle',':','LineWidth',2.5);
% vline(lthr,'color','red','LineStyle',':','LineWidth',2.5);

[pmax, pmin] = peakdet(k,peakdelta*max(k),x);

vline(pmax,'color','r','LineWidth',2);
vline(pmin,'color','g','LineWidth',2);

% [hdt_dip, hdt_p] = HartigansDipSignifTest(sprt,500);
% bc = wz_bimod_test(sprt);

xlim(spt_rng);

set(gca,'TickDir','out')
title('raw spontaneous rates')
xlabel('spike rate [spk/s]')
ylabel('count')

subplot('Position',[0.79 0.71  0.2 0.265]);
% hist(bootmean,50);
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','k','EdgeColor','k')
%
% vline(mean(bootmean),'color','b','LineWidth',2.5);
% vline(mean(bootmean) + std(bootmean),'color','b','LineStyle','--','LineWidth',2.5);
% vline(mean(bootmean) - std(bootmean),'color','b','LineStyle','--','LineWidth',2.5);

hold on;
if(sum(sprt) ~= 0)
    [t, xh] = hist(sm_sprt,25);

    bar(xh,t,1,'b','EdgeColor','b')
    [k, x] = ksdensity(sm_sprt, 'width', max(diff(xh)));
    plot(x,k*min(diff(xh))*length(sm_sprt),'-r','LineWidth',2.5);
    [pmax, pmin] = peakdet(k,peakdelta*max(k),x);

    % sm_bar = smooth(t, 20/length(xh),'loess');
    %
    % plot(xh,sm_bar,'g');
    %
    % [pmax, pmin] = peakdet(sm_bar,0.002,xh);

    vline(pmax,'color','r','LineWidth',2);
    vline(pmin,'color','g','LineWidth',2);

    [hdt_dip, hdt_p] = HartigansDipSignifTest(sm_sprt,500);
    bc = wz_bimod_test(sm_sprt);

    Ustats.sprt_dip    = hdt_dip;
    Ustats.sprt_dipP   = hdt_p;
    Ustats.sprt_bc     = bc;
    Ustats.sprt_peak   = pmax;
    Ustats.sprt_trough = pmin;

    title('filtered spontaneous rates')
    set(gca,'TickDir','out')
    cxlim = xlim;
    xlim([0,cxlim(2)]);
else
    Ustats.sprt_dip    = NaN;
    Ustats.sprt_dipP   = NaN;
    Ustats.sprt_bc     = NaN;
    Ustats.sprt_peak   = NaN;
    Ustats.sprt_trough = NaN;
end
xlabel('spike rate [spk/s]')

%% interspike
subplot('Position',[0.545 0.41  0.225 0.265]);
hold on;

if(~isempty(ISIobj.isi_hist))
    bar(ISIobj.isi_x, ISIobj.isi_hist,1,'k','EdgeColor','k');


%     if(~isempty(ISIobj.gamma_fit))
%         gf = pdf(ISIobj.gamma_fit, ISIobj.isi_x);
%         plot(ISIobj.isi_x,gf,'LineWidth',2,'color','blue')
%     end
%
%     if(~isempty(ISIobj.lognorm_fit))
%         lf = pdf(ISIobj.lognorm_fit, ISIobj.isi_x);
%         plot(ISIobj.isi_x,lf,'LineWidth',2,'color','red')
%     end

    vline(mean(ISIobj.isi_vec),  'color','green');
    vline(median(ISIobj.isi_vec),'color','cyan');

    xlim([0,prctile(ISIobj.isi_vec, 75)]);
    set(gca,'TickDir','out')
end
xlabel('interspike interval [ms]');

subplot('Position',[0.79 0.41 0.2 0.265]);
hold on;

if(~isempty(ISIobj.isi_hist))
[t, xh] = hist(sm_isi,25);

    bar(xh,t,1,'b','EdgeColor','b')
    [k, x] = ksdensity(sm_isi, 'width', max(diff(xh)));
    plot(x,k*min(diff(xh))*length(sm_sprt),'-r','LineWidth',2.5);

    [pmax, pmin] = peakdet(k,peakdelta*max(k),x);

    vline(pmax,'color','r','LineWidth',2);
    vline(pmin,'color','g','LineWidth',2);

    [hdt_dip, hdt_p] = HartigansDipSignifTest(sm_isi,500);
    bc = wz_bimod_test(sm_isi);

    Ustats.isi_dip    = hdt_dip;
    Ustats.isi_dipP   = hdt_p;
    Ustats.isi_bc     = bc;
    Ustats.isi_peak   = pmax;
    Ustats.isi_trough = pmin;

    % title('filtered ISI')
    set(gca,'TickDir','out')
    xlabel('spike rate [spk/s]')
    cxlim = xlim;
    if(cxlim(2) > 0)
        xlim([0,cxlim(2)]);
    end
else
    Ustats.isi_dip    = NaN;
    Ustats.isi_dipP   = NaN;
    Ustats.isi_bc     = NaN;
    Ustats.isi_peak   = NaN;
    Ustats.isi_trough = NaN;
end

xlabel('interspike interval [ms]')

subplot('Position',[0.545 0.075  0.325 0.25]);
hold on;

if(sum(isfinite(spontspk(:))) > 10)
    [ov, xv] = wz_spk_RACH(spontspk, [], 150);

    bar(xv, ov,1,'k','EdgeColor','k');
    % [k, x] = ksdensity(acval, 'bandwidth', 2);
    % plot(x,k*length(acval),'-r','LineWidth',2);

    % sm_ac = smooth(ov, 20/length(xv),'loess');
    % plot(xv, sm_ac,'-r','LineWidth',2);

    if(any(ov ~= 0))
        acv = sort(ov);
        ylim([0.9*acv(2),1.05*acv(end)]);
    end

    xlim([-150, 150]);
    xlabel('lag [ms]')
    set(gca,'TickDir','out')
end

%% find outliers in spontaneous rates
uthr = median(sm_sprt) + 2*mad(sm_sprt);
lthr = median(sm_sprt) - 2*mad(sm_sprt);

% identify trials with very different spontaneous rates
oltrials = find(sm_sprt > uthr | sm_sprt < lthr);
blkpos   = strfind(diff(oltrials)', ones(1,olblk-1));

blckend = 0;
outblock = [];
for(i=1:length(blkpos))
    p = find(diff(oltrials(blkpos(i):end))~=1,1,'first');
    if(p<blckend)
        continue;
    end
    if(isempty(p))
        blckend = spk.nTrials;
    else
        blckend = oltrials(blkpos(i)+p-1);
    end
    outblock = [outblock; oltrials(blkpos(i)), blckend];
end

% sp_oltrials           = oltrials;
Ustats.sprt_outblock  = outblock;

uthr = median(sm_isi) + 2*mad(sm_isi);
lthr = median(sm_isi) - 2*mad(sm_isi);

% identify trials with very different spontaneous rates
oltrials = find(sm_isi > uthr | sm_isi < lthr);
blkpos   = strfind(diff(oltrials)', ones(1,olblk-1));

blckend = 0;
outblock = [];
for(i=1:length(blkpos))
    p = find(diff(oltrials(blkpos(i):end))~=1,1,'first');
    if(p<blckend)
        continue;
    end
    if(isempty(p))
        blckend = spk.nTrials;
    else
        blckend = oltrials(blkpos(i)+p-1);
    end
    outblock = [outblock; oltrials(blkpos(i)), blckend];
end

% isi_oltrials        = oltrials;
Ustats.isi_outblock = outblock;

%% look for steps of firing rate
Ustats.rate_steps = [];
Ustats.rate_p     = [];
i=minbblock;

while(i<spk.nTrials-minbblock)
    preblk = i-minbblock+1:i-1;
    medrt  = median(sm_sprt(preblk));
    madrt  = (median(abs(sm_sprt(preblk) - medrt)));

    if(~any(abs(sm_sprt(i:i+round(minbblock/2))-medrt) < 2*madrt)) % transition possible
        rtP = ranksum(sprt(preblk),sprt(i:i+minbblock-1));
        if(rtP <= 0.001)    % transition
            Ustats.rate_steps = [Ustats.rate_steps , i];
            Ustats.rate_p     = [Ustats.rate_p, rtP];
            i = i + minbblock;
        else
            i = i + 1;
        end
    else
        i = i + 1;
    end
end

Ustats.rate_Ngroup = length(Ustats.rate_steps)+1;
Ustats.rate_groups = ones(1,spk.nTrials);

for(i=1:length(Ustats.rate_steps))
    Ustats.rate_groups(Ustats.rate_steps(i):end) = i+1;
end

if(Ustats.rate_Ngroup> 1)
    rep = 5000;
    smplsz = 25;

    rate_fdist   = nan(rep,1);
    rate_pdist   = nan(rep,1);
    rate_smplmat = [];

    grpvec  = [];
    for(i=1:Ustats.rate_Ngroup)
        gpos = find(Ustats.rate_groups == i);
        rate_smplmat = [rate_smplmat;  sprt(gpos(randi(length(gpos), smplsz, rep)))];
        grpvec  = [grpvec; repmat(i,smplsz,1)];
    end

    for(r=1:rep)
        [~,aovtbl]    = anova1(rate_smplmat(:,r), grpvec,'off');
        rate_fdist(r) = cell2mat(aovtbl(2,5));
        rate_pdist(r) = cell2mat(aovtbl(2,6));
    end

    Ustats.group_rate_meanP   = mean(  rate_pdist);
    Ustats.group_rate_meanF   = mean(  rate_fdist);
    Ustats.group_rate_medianP = median(rate_pdist);
    Ustats.group_rate_medianF = median(rate_fdist);

else
    Ustats.group_rate_meanP   = NaN;
    Ustats.group_rate_meanF   = NaN;
    Ustats.group_rate_medianP = NaN;
    Ustats.group_rate_medianF = NaN;
end

%% determine stepchange
isi_grp = ones(spk.nTrials,1);
if(~isempty(Ustats.isi_trough))
    cgrp=2;
    for(i=1:length(Ustats.isi_trough)-1)
        p  = sm_isi > Ustats.isi_trough(i) & sm_isi > Ustats.isi_trough(i+1);
        isi_grp(p) = cgrp;
        cgrp = cgrp + 1;
    end
    isi_grp(sm_isi > Ustats.isi_trough(end)) = cgrp;
end
Ustats.isi_grp = isi_grp;

sprt_grp = ones(spk.nTrials,1);
if(~isempty(Ustats.sprt_trough))
    cgrp=2;
    for(i=1:length(Ustats.sprt_trough)-1)
        p  = sm_sprt > Ustats.sprt_trough(i) & sm_sprt > Ustats.sprt_trough(i+1);
        sprt_grp(p) = cgrp;
        cgrp = cgrp + 1;
    end
    sprt_grp(sm_sprt > Ustats.sprt_trough(end)) = cgrp;
end
Ustats.sprt_grp = sprt_grp;

%% use cluster to partion response states
if(~isempty(Ustats.sprt_trough) || ~isempty(Ustats.isi_trough) || Ustats.sprt_bc > 0.5 || Ustats.isi_bc > 0.5 || Ustats.isi_dipP < 0.05 ||  Ustats.sprt_dipP < 0.05)

    ntrough = max([length(Ustats.sprt_trough), length(Ustats.isi_trough)]);
    if(ntrough+1>2)
        ngrps = ntrough+1;
    else
        ngrps = 2;
    end

    trn = 1:spk.nTrials;

    mad_sm_spk = (sm_sprt - median(sm_sprt)) ./ mad(sm_sprt);
    mad_sm_isi = (sm_isi  - median(sm_isi))  ./ mad(sm_isi);

%     distvec = sqrt(mad_sm_spk.^2 +  mad_sm_isi.^2);
   tm = [sm_sprt(:), sm_isi(:), trn(:)];
%    tm = [mad_sm_spk(:), mad_sm_isi(:), trn(:)];

    try
        obj = gmdistribution.fit(tm,ngrps);
        [idx,nlogl] = cluster(obj, tm);
        idx = idx';
    catch
        idx = ones(1,spk.nTrials);
        nlogl = NaN;
    end

    statetrans = find(diff(idx)~=0);

    statetrans(statetrans<50) = [];

    if(~isempty(statetrans))
    % identify blocks that are too short to be accounted for
        p = find(diff(statetrans) < minbblock);

        while(~isempty(p))
            if(~isempty(p))
                for(i=length(p):-1:1)
                   if(idx(statetrans(p(i))-1) == idx(statetrans(p(i)+1)+1))
                        idx(statetrans(p(i)):statetrans(p(i)+1)+1) = idx(statetrans(p(i))-1);
                   else
                       pBef = idx == idx(statetrans(p(i))-1);
                       pAft = idx == idx(statetrans(p(i)+1)+1);
                       pSmp = idx == idx(statetrans(p(i))+1);

                       cdistB = sqrt((mean(mad_sm_spk(pBef))-mean(mad_sm_spk(pSmp)))^2 ...
                                   + (mean(mad_sm_isi(pBef))-mean(mad_sm_isi(pSmp)))^2);
                       cdistA = sqrt((mean(mad_sm_spk(pAft))-mean(mad_sm_spk(pSmp)))^2 ...
                                   + (mean(mad_sm_isi(pAft))-mean(mad_sm_isi(pSmp)))^2);

                       if(cdistA < cdistB)
                           idx(statetrans(p(i)):statetrans(p(i)+1)+1) = idx(statetrans(p(i)+1)+1);
                       else
                           idx(statetrans(p(i)):statetrans(p(i)+1)+1) = idx(statetrans(p(i))-1);
                       end
                   end
                end

                if(statetrans(1) < minbblock)
                    idx(1:statetrans(1)+1) = idx(statetrans(1)+1);
                end
                if(length(idx) - statetrans(end) < minbblock)
                    idx(statetrans(end):end) = idx(statetrans(end)-1);
                end

                idxlst = unique(idx);
                nidx = nan(1,length(idx));
                for(i=1:length(idxlst))
                    nidx(idx==idxlst(i)) = i;
                end
                idx = nidx;
            end
            statetrans = find(diff(idx)~=0);  % check again
            statetrans(statetrans<50) = [];
            p = find(diff(statetrans) < minbblock);
        end  %  while(~isempty(p))
    end  %  if(~isempty(statetrans))

    % order block indices
    if(~isempty(statetrans))
        nidx = nan(1,length(idx));
        sblock = 1;
        cblock = 1;
        for(i=1:length(statetrans))
            nidx(sblock:statetrans(i)) = cblock;
            sblock = statetrans(i) + 1;
            cblock = cblock + 1;
        end
        nidx(sblock:end) = cblock;
        idx = nidx;
    end

    Ustats.trial_groups       = idx;
    Ustats.trial_groups_nlogl = nlogl;

    Ustats.trial_group_trans = statetrans;
    Ustats.trial_Ngroup      = length(unique(idx));

    for(i=1:Ustats.trial_Ngroup)
        Ustats.trial_groupSZ(i) = sum(idx == i);
    end
    Ustats.trial_groupProp = Ustats.trial_groupSZ ./ length(idx);

else
    Ustats.trial_groups       = ones(spk.nTrials,1);
    Ustats.trial_groups_nlogl = NaN;
    Ustats.trial_group_trans  = [];
    Ustats.trial_Ngroup       = 1;
    Ustats.trial_groupSZ      = spk.nTrials;
    Ustats.trial_groupProp    = 1;
end

if(Ustats.trial_Ngroup > 1)
    rep = 5000;
    smplsz = 25;

    spk_fdist = nan(rep,1);
    spk_pdist = nan(rep,1);
    spk_smplmat = [];

    isi_fdist = nan(rep,1);
    isi_pdist = nan(rep,1);
    isi_smplmat = [];

    grpvec  = [];
    for(i=1:Ustats.trial_Ngroup)
        gpos = find(Ustats.trial_groups == i);
        spk_smplmat = [spk_smplmat;  sprt(gpos(randi(length(gpos), smplsz, rep)))];
        isi_smplmat = [isi_smplmat;  ISIobj.trial_isi(gpos(randi(length(gpos), smplsz, rep)))];
        grpvec  = [grpvec; repmat(i,smplsz,1)];
    end

    for(r=1:rep)
        [~,aovtbl]   = anova1(spk_smplmat(:,r), grpvec,'off');
        spk_fdist(r) = cell2mat(aovtbl(2,5));
        spk_pdist(r) = cell2mat(aovtbl(2,6));

        [~,aovtbl]   = anova1(isi_smplmat(:,r), grpvec,'off');
        isi_fdist(r) = cell2mat(aovtbl(2,5));
        isi_pdist(r) = cell2mat(aovtbl(2,6));
    end

    Ustats.group_spk_meanP   = mean(  spk_pdist);
    Ustats.group_spk_meanF   = mean(  spk_fdist);
    Ustats.group_spk_medianP = median(spk_pdist);
    Ustats.group_spk_medianF = median(spk_fdist);

    Ustats.group_isi_meanP   = mean(  isi_pdist);
    Ustats.group_isi_meanF   = mean(  isi_fdist);
    Ustats.group_isi_medianP = median(isi_pdist);
    Ustats.group_isi_medianF = median(isi_fdist);
else
    Ustats.group_spk_meanP   = NaN;
    Ustats.group_spk_meanF   = NaN;
    Ustats.group_spk_medianP = NaN;
    Ustats.group_spk_medianF = NaN;

    Ustats.group_isi_meanP   = NaN;
    Ustats.group_isi_meanF   = NaN;
    Ustats.group_isi_medianP = NaN;
    Ustats.group_isi_medianF = NaN;
end

%% plot rate time course
subplot('Position',[0.305 0.075  0.1 0.9]);
hold on;

trial_vec = 1:spk.nTrials;

grpcol = [155, 255, 255; 255, 255, 155; 255, 155, 255; ...
          155, 155, 200; 155, 200, 155; 200, 155, 155; ...
          230, 255, 200; 230, 200, 255; 255, 230, 200; 255, 200, 230; 200, 230, 255; 200, 255, 230]./255;
lend = 1;
for(i=1:length(Ustats.rate_steps))
    lstart = Ustats.rate_steps(i);
    patch(spt_rng([1 1 2 2]), [lstart lend lend lstart], grpcol(mod(i,12)+1,:), 'EdgeColor', grpcol(mod(i,12)+1,:));
    lend = lstart+1;
end

% for(i=1:size(outblock,1))
%     patch(spt_rng([1 1 2 2]), outblock(i,[1 2 2 1]), 'yellow', 'EdgeColor', 'yellow');
% end

plot(sprt, trial_vec,'k','LineWidth',0.8);

%plot(sprt(sp_oltrials), trial_vec(sp_oltrials),'.m');

plot(sm_sprt, trial_vec, 'b', 'LineWidth', 2.5);

% vline(spt_med,'color','r','LineWidth',2.5);
% vline(spt_med+spt_mad,'color','r','LineWidth',2,'LineStyle','--');
% vline(spt_med-spt_mad,'color','r','LineWidth',2,'LineStyle','--');
% vline(uthr,'color','r','LineWidth',2,'LineStyle',':');
% vline(lthr,'color','r','LineWidth',2,'LineStyle',':');

vline(median(sm_sprt) + 2*mad(sm_sprt,1),'color','b','LineWidth',2,'LineStyle','--');
vline(median(sm_sprt) - 2*mad(sm_sprt,1),'color','b','LineWidth',2,'LineStyle','--');

hline(Ustats.rate_steps,'color','r','LineWidth',2);

set(gca,'YTick',[]);
set(gca,'TickDir','out')
ylim([0; spk.nTrials+1]);
xlim(spt_rng);
xlabel('spike rate [spk/s]')

subplot('Position',[0.41 0.075  0.1 0.9]);
hold on;

isi_rng = [min(ISIobj.trial_isi), max(ISIobj.trial_isi)];

lend = 1;
for(i=1:length(Ustats.trial_group_trans))
    lstart = Ustats.trial_group_trans(i);
    patch(isi_rng([1 1 2 2]), [lstart lend lend lstart], grpcol(mod(i,12)+1,:), 'EdgeColor', grpcol(mod(i,12)+1,:));
    lend = lstart+1;
end

plot(ISIobj.trial_isi, trial_vec,'k','LineWidth',0.8);

%plot(sprt(sp_oltrials), trial_vec(sp_oltrials),'.m');

plot(sm_isi, trial_vec, 'b', 'LineWidth', 2.5);

vline(median(sm_isi) + 2*mad(sm_isi,1),'color','b','LineWidth',2,'LineStyle','--');
vline(median(sm_isi) - 2*mad(sm_isi,1),'color','b','LineWidth',2,'LineStyle','--');

hline(Ustats.trial_group_trans,'color','c','LineWidth',2);

set(gca,'YTick',[]);
set(gca,'TickDir','out')
ylim([0; spk.nTrials+1]);
xlim(isi_rng);
xlabel('mean ISI [ms]')


%% write data
subplot('Position',[0.89 0.075  0.115 0.3]);
hold on;
axis off

smry_str = {[nfo], ...
            ['mean rate = ', num2str(mean(sprt),'%.2f'), 'spks/s'] ...
            ['mean ISI = ', num2str(mean(sprt),'%.2f'), 'ms'] ...
            [''],...
            ['First Bins / Peak = ', num2str(Ustats.isi1stBins2Peak,'%.2f')], ...
            ['First Bins / All = ', num2str(Ustats.isi1stBins2all,'%.2f')], ...
            ['ISI half Max Bin = ', int2str(Ustats.ISIHalfMaxBin)], ...
            [''],...
            ['Spont Rate BC = ', num2str(Ustats.sprt_bc,'%.4f')], ...
            ['Spont Rate dip = ', num2str(Ustats.sprt_dip,'%.4f')], ...
            ['Spont Rate dip P = ', num2str(Ustats.sprt_dipP,'%.4f')], ...
            ['ISI BC = ', num2str(Ustats.isi_bc,'%.4f')], ...
            ['ISI dip = ', num2str(Ustats.isi_dip,'%.4f')], ...
            ['ISI dip P = ', num2str(Ustats.isi_dipP,'%.4f')], ...
            [''],...
            ['N Blocks = ', int2str(Ustats.trial_Ngroup)], ...
            ['Blocks nlogl = ', num2str(Ustats.trial_groups_nlogl,'%.2f ')], ...
            ['Block Size = ', int2str(Ustats.trial_groupSZ)], ...
            ['Block Proportion = ', num2str(Ustats.trial_groupProp,'%.2f ')], ...
            ['Block Transitions = ', int2str(Ustats.trial_group_trans)], ...
            [''],...
            ['Rate Transitions = ', int2str(Ustats.rate_steps)], ...
            ['Rate P = ', num2str(Ustats.rate_p,'%.2f ')], ...
            ['Rate median F = ', num2str(Ustats.group_rate_medianF,'%.4f')], ...
            ['Rate median P = ', num2str(Ustats.group_rate_medianP,'%.4f')], ...
            [''],...
            ['Spont median F = ', num2str(Ustats.group_spk_medianF,'%.4f')], ...
            ['Spont median P = ', num2str(Ustats.group_spk_medianP,'%.4f')], ...
            ['ISI median F = ', num2str(Ustats.group_isi_medianF,'%.4f')], ...
            ['ISI median P = ', num2str(Ustats.group_isi_medianP,'%.4f')]};

text(0, 1, 1, smry_str, 'FontSize', 10, 'VerticalAlignment', 'top');
