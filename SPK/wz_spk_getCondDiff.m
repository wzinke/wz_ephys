function TST = wz_spk_getCondDiff(spks, MGRF, condvec, cliptime, timwin, spontwin, nBoot, smplsz, krnl, kw, doplot, nfo)
% ToDo: be less conservative with the Time window limit for comparisions of two conditions.

% %% check input and define default values
% if(~isfield(spk1,'spikedensities'))
%     spk1 = wz_spk_density(spk1, 'exp', 20);
% end

% p = dt.correct == 1 & dt.SetSZ == 4;
% %p = dt.correct == 1 & dt.SetSZ == 8;
% condvec = dt.Tpos(p);
% spks = dt.spiketimes(p,:);
% cliptime = dt.SRT(p);
%
% MGRF = dt.RF;

%  ========================================================================
%% Define constant parameter

theta_vals = [0:45:315];
% thv        = theta_vals ./ (180/pi);
% angstep    = 45/(180/pi);

mincons = 20;     % at least <mincons> time points below the <siglev> to be considered

%  ========================================================================
%% check input arguments and define defaults

% number of bootstrap repetitions
if(~exist('nBoot','var') || isempty(nBoot))
    nBoot = 1;
end

% number of bootstrap repetitions
if(~exist('smplsz','var') || isempty(smplsz))
    if(nBoot == 1)
        smplsz = [];
    else
        smplsz = 25;
    end
end

if(~exist('timwin','var') || isempty(timwin))
    timwin = [30, 300];
end

if(~exist('spontwin','var') || isempty(spontwin))
    spontwin = [-500, 0];
end

if(~exist('MGRF','var'))
    MGRF = [];
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

% =========================================================================
%% prepare data
% spksclip = spks(bsxfun(@lt, spks, cliptime));
% spksclip = wz_cropNaN(spksclip, [], timwin(2));
spksclip = spks;
spksclip(spksclip > timwin(2)) = NaN;
spksclip(bsxfun(@gt, spksclip, cliptime)) = NaN;

respraw = wz_spk_density(spks,     krnl, kw);
respest = wz_spk_density(spksclip, krnl, kw, 1, cliptime); % remove spikes that occurred after responses

% identify time bins with meaningful values
respest.spikedensities(bsxfun(@gt, respest.xtime, cliptime)) = NaN;
trialmat = isfinite(respest.spikedensities);

% get list of available conditions
cndvec = unique(condvec(condvec ~= 255));
cndcnt = hist(condvec, cndvec); % how many trials per condition?

if(isempty(smplsz))
%     smplsz = min(cndcnt);
    smplsz = round(median(cndcnt)); % not so clean ad hoc fix to cope with situations where one condition has no data.
end

if(isempty(smplsz))
    smplsz1 = Ntrials_1;
    smplsz2 = Ntrials_2;
else
    smplsz1 = smplsz;
    smplsz2 = smplsz;
end

% get index of relevant time bins
xpos = respest.xtime >= timwin(1) & respest.xtime < timwin(2);
xvec = respest.xtime(xpos);
xrng = [min(xvec), max(xvec)];

% ----------------------------------------------------------------------- %
% get mean response estimate
meanraw     = nan(length(cndvec), size(respraw.spikedensities,2));
errraw      = nan(length(cndvec), size(respraw.spikedensities,2));

meanresp    = nan(length(cndvec), size(respest.spikedensities,2));
errresp     = nan(length(cndvec), size(respest.spikedensities,2));
rateresp    = nan(length(cndvec), 1);
cndSRTquart = nan(length(cndvec), 3);

for(i=1:length(cndvec))
    cpos = condvec == cndvec(i);
    meanraw(i,:)   = nanmean(   respraw.spikedensities(cpos, :));
    trialcnt       = sum(~isnan(respraw.spikedensities(cpos, :)));
    errraw(i,:)    = nanstd(    respraw.spikedensities(cpos, :),1) ./ sqrt(trialcnt);

    meanresp(i,:)  = nanmean(respest.spikedensities(cpos, :));
    trialcnt       = sum(trialmat(cpos, :));
    errresp(i,:)   = nanstd(respest.spikedensities(cpos, :),1) ./ sqrt(trialcnt);

    rateresp(i)    = nanmean(nansum(respest.spikedensities(cpos, xpos),2) ./ cliptime(cpos));

    cndSRTquart(i,:) = prctile(cliptime(cpos),[25 50 75]);
end

% ----------------------------------------------------------------------- %
% identify the preferred locations and its counterpart
% respmax = max(meanresp,[],2);
% respmin = min(meanresp,[],2);

% determine locations with driven responses
spontpos  = respest.xtime >= spontwin(1) & respest.xtime < spontwin(2);
% spontrate = 1000 * (nansum(respraw.spikedensities(:, spontpos),2) ./ diff(spontwin));

% VisResp = find(respmax > nanmedian(spontrate) + 2*mad(spontrate));
% InRF    = unique(mod(VisResp+3,8)+1);

prefloc = find(rateresp == max(rateresp));
antiloc = mod(prefloc+3,8)+1;

prefPos = find(condvec == mod(prefloc,8));
antiPos = find(condvec == mod(antiloc,8));

if(~isempty(MGRF) && length(MGRF) <= 4)
    InRFpos  = find(ismember(condvec, MGRF) == 1);
%     OutRFpos = find(ismember(condvec, mod(MGRF+3,8)+1) == 1);

    OutRFpos = find(ismember(condvec, mod(MGRF+4,8)) == 1); % zero based? DOUBLECHECK!!!

else
    InRFpos  = [];
    OutRFpos = [];
end

% prefRad = circ_mean(thv', rateresp);
% prefDir = prefRad / (pi/180);

% ========================================================================
%% create matrix for bootstrap repetitions
rng(sum(100*clock), 'twister');

if(nBoot>1)
    Bsmpl = nan(nBoot,smplsz*length(cndvec));
    for(i=1:length(cndvec))
        cpos = find(condvec == cndvec(i));
        cs = (i-1)*smplsz+1;
        Bsmpl(:,cs:cs+smplsz-1) = cpos(randi(length(cpos),  nBoot, smplsz));
    end

    if(~isempty(prefPos) && ~isempty(antiPos))
        MXSmplT = randi(length(prefPos), nBoot, smplsz);
        MXSmplD = randi(length(antiPos), nBoot, smplsz);
    end

    if(~isempty(InRFpos) && ~isempty(OutRFpos) && length(MGRF) <= 4)
        RFsmplT  = randi(length(InRFpos),  nBoot, smplsz);
        RFsmplD  = randi(length(OutRFpos), nBoot, smplsz);
    end
else
    Bsmpl = 1:length(condvec);

    MXSmplT = 1:length(prefPos);
    MXSmplD = 1:length(antiPos);
    RFsmplT = 1:length(InRFpos);
    RFsmplD = 1:length(OutRFpos);
end

% ========================================================================
%% get a bootstrap estimate of response times
% bSRT_all = nan(nBoot, s1);
%
%
% bSRT_MX    = nan(nBoot, s1);
% bSRT_MXin  = nan(nBoot, s1);
% bSRT_MXout = nan(nBoot, s1);
%
% for(b=1:nBoot) % loop over bootstrap repetitions
%     csrt = cliptime(Bsmpl(b,:));
%     cloc = condvec( Bsmpl(b,:));
%
%     bSRT_all(b,1) = mean(csrt);
%
%     bSRT_RF    = nan(nBoot, s1);
%
%     InRFpos
%
%     cTres = inRF_T(RFsmplT(b,:));
%     cTres(isnan(cTres)) = [];
%     cDres = inRF_D(RFsmplD(b,:));
%     cDres(isnan(cDres)) = [];
%
%     cTres = MX_T(MXSmplT(b,:));
%     cTres(isnan(cTres)) = [];
%     cDres = MX_D(MXSmplD(b,:));
%     cDres(isnan(cDres)) = [];
% end

% ========================================================================
%% run the estimations by testing all conditions
pvec    = nan(nBoot, sum(xpos));
% bSRT_RF = nan(nBoot, 1);

% ensure that there is sufficient data for each location
if(timwin(2) > min(cndSRTquart(:,3)))
    numtestbins = ceil(min(cndSRTquart(:,3)) - timwin(1));
    if(numtestbins > sum(xpos))
       numtestbins = sum(xpos);
    end
else
    numtestbins = sum(xpos);
end

spkdense = respest.spikedensities(:, xpos);

%     bSRT_all(b,1) = mean(cliptime(Bsmpl(b,:)));

parfor(i=1:numtestbins) % loop over time bins
    for(b=1:nBoot) % loop over bootstrap repetitions
        csr  = spkdense(Bsmpl(b,:),i);
        cloc = condvec( Bsmpl(b,:));

        cloc(isnan(csr)) = [];
        csr( isnan(csr)) = [];

        % [pval, table, stats] = anova1(csr,cloc,'off');
        pvec(b,i) = kruskalwallis(csr, cloc, 'off');
    end
end

% ----------------------------------------------------------------------- %
% get estimate for target within RF and outside RF
pRF = nan(nBoot, sum(xpos));

% bSRT_RFall = nan(nBoot, 1);
% bSRT_RFin  = nan(nBoot, 1);
% bSRT_RFout = nan(nBoot, 1);

if(~isempty(InRFpos) && ~isempty(OutRFpos) && length(MGRF) <= 4)
    Tcnt = size(RFsmplT,1);  % use this measure to ensure that at least 75% valid values are used
    Dcnt = size(RFsmplD,1);

%         bSRT_RFall(b,1) = mean(cliptime(InRFpos(b,:) + ));
%
%         bSRT_RFin  = mean(cliptime(InRFpos( RFsmplT(b,:))));
%         bSRT_RFout = mean(cliptime(OutRFpos(RFsmplD(b,:))));

    parfor(i=1:sum(xpos)) % loop over time bins
        inRF_T = spkdense(InRFpos, i);
        inRF_D = spkdense(OutRFpos, i);

        for(b=1:nBoot) % loop over bootstrap repetitions

            cTres = inRF_T(RFsmplT(b,:));
            cTres(isnan(cTres)) = [];
            cDres = inRF_D(RFsmplD(b,:));
            cDres(isnan(cDres)) = [];

            if(length(cTres) >= Tcnt/4 && length(cDres) >= Dcnt/4)
                pRF(b,i) = ranksum(cTres, cDres, 'method','approximate');
            end
        end
    end
end

% ----------------------------------------------------------------------- %
% get estimate for difference of preferred and opposite location
pMax = nan(nBoot, sum(xpos));

if(~isempty(prefloc))
    if(~isempty(prefPos) && ~isempty(antiPos))
        Tcnt = size(MXSmplT,1);  % use this measure to ensure that at least 75% valid values are used
        Dcnt = size(MXSmplD,1);

        parfor(i=1:sum(xpos)) % loop over time bins
            MX_T = spkdense(prefPos, i);
            MX_D = spkdense(antiPos, i);

            for(b=1:nBoot) % loop over bootstrap repetitions
                cTres = MX_T(MXSmplT(b,:));
                cTres(isnan(cTres)) = [];
                cDres = MX_D(MXSmplD(b,:));
                cDres(isnan(cDres)) = [];

                if(length(cTres) >= Tcnt/4 && length(cDres) >= Dcnt/4)
                    pMax(b,i) = ranksum(cTres, cDres, 'method','approximate');
                end
            end
        end
    end
end

% ----------------------------------------------------------------------- %
% utilize circular statistics to derive tuning properties
% tv   = nan(nBoot, length(xpos));
% tdir = nan(nBoot, length(xpos));
% pur  = nan(nBoot, length(xpos));
% puo  = nan(nBoot, length(xpos));
% puv  = nan(nBoot, length(xpos));
% angvec = thv(condvec+1);
% for(i=1:length(xpos)) % loop over time bins
%     tp = xpos(i);
%     bin_sr = respest.spikedensities(:, tp);
%
%     parfor(b=1:nBoot) % loop over bootstrap repetitions
%         csr  = bin_sr( Bsmpl(b,:));
%         cloc = condvec(Bsmpl(b,:));
%
%         cloc(isnan(csr)) = [];
%         csr( isnan(csr)) = [];
%
% %         pvec(b,i) = pval;
%
% %         angvec = thv(cloc+1);
% %         tv(b,i)   = circ_r(angvec', csr./max(csr), angstep);
% %         tdir(b,i) = circ_mean(angvec', csr);
%     end
% end

%  ========================================================================
%% get estimate based on response differences
DiffTDTpref     = NaN;
bootDiffTDTpref = nan(1,nBoot);
if(~isempty(prefPos) && ~isempty(antiPos))

    [PrefBdiff, PrefBstd, bootdiff] = SPK_boot_RespDiff(respest.spikedensities(prefPos,:), ...
                                                        respest.spikedensities(antiPos,:), nBoot, smplsz);
    PrefBdiff = abs(PrefBdiff);
    DiffXced  = PrefBdiff(xpos) > prctile(PrefBdiff(spontpos)+PrefBstd(spontpos), 99);
    sigpos    = min(strfind(DiffXced, ones(1,mincons)));

    if(~isempty(sigpos))
        DiffTDTpref  = xvec(sigpos);
    else
        DiffTDTpref  = NaN;
    end

    bootdiff = abs(bootdiff);

    DiffXced_mat = bsxfun(@gt, bootdiff(:,xpos), prctile(bootdiff(:,spontpos), 99,2));

    parfor(r=1:nBoot)
        sigpos = min(strfind(DiffXced_mat(r,:), ones(1,mincons)));
        if(~isempty(sigpos))
            bootDiffTDTpref(r) = xvec(sigpos);
        end
    end
end

%  ========================================================================
DiffTDTrf = NaN;
MaxBdiff  = NaN;
MaxBstd   = NaN;

bootDiffTDTrf = nan(1,nBoot);

if(~isempty(InRFpos) && ~isempty(OutRFpos) && length(MGRF) <= 4)
    [MaxBdiff, MaxBstd, bootdiff] = SPK_boot_RespDiff(respest.spikedensities(InRFpos, :), ...
                                                      respest.spikedensities(OutRFpos, :), nBoot, smplsz);
    MaxBdiff = abs(MaxBdiff);
    DiffXced = MaxBdiff(xpos) > prctile(MaxBdiff(spontpos)+MaxBstd(spontpos), 99);
    sigpos   = min(strfind(DiffXced, ones(1,mincons)));
    if(~isempty(sigpos))
        DiffTDTrf  = xvec(sigpos);
    end

    bootdiff = abs(bootdiff);
    DiffXced_mat = bsxfun(@gt, bootdiff(:,xpos), prctile(bootdiff(:,spontpos), 99,2));

    parfor(r=1:nBoot)
        sigpos = min(strfind(DiffXced_mat(r,:), ones(1,mincons)));
        if(~isempty(sigpos))
            bootDiffTDTrf(r) = xvec(sigpos);
        end
    end
end

%  ========================================================================
%% get tuning and p time courses
median_pvec = nanmedian(pvec);
upperP      = prctile(  pvec,75,1);
lowerP      = prctile(  pvec,25,1);

median_pMax = nanmedian(pMax);
upperPMax   = prctile(  pMax,75,1);
lowerPMax   = prctile(  pMax,25,1);

median_pRF  = nanmedian(pRF);
upperPRF    = prctile(  pRF,75,1);
lowerPRF    = prctile(  pRF,25,1);

%mean_tune = nanmean(tv);
% std_tune  =  nanstd(tv);
%
% mean_dir    = circ_mean(tdir);
% [ad, cs]    = circ_std( tdir);
% AngMean     = rad2deg(abs(mean_dir));
% AngDev      = rad2deg(ad);
% AngStd      = rad2deg(cs);
% p = isfinite(mean_tune) & isfinite(mean_dir);
% mndir        = circ_mean(mean_dir(p), mean_tune(p), 2);
% [mnad, mncs] = circ_std( mean_dir(p), mean_tune(p), [], 2);
% mnAng      = rad2deg(abs(mndir));
% mnAngDev   = rad2deg(mnad);
% mnAngStd   = rad2deg(mncs);

%  ========================================================================
%%  determine time point of clear response difference
[TDT05,  TDT05conf]  = TDT_getsigdiff(median_pvec, upperP, lowerP, xvec, 0.05,  mincons);
[TDT01,  TDT01conf]  = TDT_getsigdiff(median_pvec, upperP, lowerP, xvec, 0.01,  mincons);
[TDT001, TDT001conf] = TDT_getsigdiff(median_pvec, upperP, lowerP, xvec, 0.001, mincons);

[TDTmax05,  TDTmax05conf]  = TDT_getsigdiff(median_pMax, upperPMax, lowerPMax, xvec, 0.05,  mincons);
[TDTmax01,  TDTmax01conf]  = TDT_getsigdiff(median_pMax, upperPMax, lowerPMax, xvec, 0.01,  mincons);
[TDTmax001, TDTmax001conf] = TDT_getsigdiff(median_pMax, upperPMax, lowerPMax, xvec, 0.001, mincons);

[TDTrf05,   TDTrf05conf]  = TDT_getsigdiff(median_pRF, upperPRF, lowerPRF, xvec, 0.05,  mincons);
[TDTrf01,   TDTrf01conf]  = TDT_getsigdiff(median_pRF, upperPRF, lowerPRF, xvec, 0.01,  mincons);
[TDTrf001,  TDTrf001conf] = TDT_getsigdiff(median_pRF, upperPRF, lowerPRF, xvec, 0.001, mincons);

%  ========================================================================
% do it on the bootstrapped data
bootTDT05     = TDT_getsigP(pvec, xvec, 0.05, mincons);
bootTDT01     = TDT_getsigP(pvec, xvec, 0.01, mincons);
bootTDT001    = TDT_getsigP(pvec, xvec, 0.001, mincons);

bootTDTmax05  = TDT_getsigP(pMax, xvec, 0.05, mincons);
bootTDTmax01  = TDT_getsigP(pMax, xvec, 0.01, mincons);
bootTDTmax001 = TDT_getsigP(pMax, xvec, 0.001, mincons);

bootTDTrf05   = TDT_getsigP(pRF, xvec, 0.05, mincons);
bootTDTrf01   = TDT_getsigP(pRF, xvec, 0.01, mincons);
bootTDTrf001  = TDT_getsigP(pRF, xvec, 0.001, mincons);

%  ========================================================================
%% define output
if(nargout>0)
    TST.nBoot   = nBoot;
    TST.smplsz  = smplsz;

    TST.prefloc = prefloc-1;  % zero based
    TST.antiloc = antiloc-1;

    TST.TDT05   = TDT05;
    TST.TDT05ci = TDT05conf;

    % response difference based
    TST.DiffTDTrf        = DiffTDTrf;
    TST.bootDiffTDTrf    = median(bootDiffTDTrf(isfinite(bootDiffTDTrf)));
    TST.bootDiffTDTrfiqr =    iqr(bootDiffTDTrf(isfinite(bootDiffTDTrf)));
    TST.bootDiffTDTrfacc = 100*sum(isfinite(bootDiffTDTrf))/length(bootDiffTDTrf);

    TST.DiffTDTpref        = DiffTDTpref;
    TST.bootDiffTDTpref    = median(bootDiffTDTpref(isfinite(bootDiffTDTpref)));
    TST.bootDiffTDTprefiqr =    iqr(bootDiffTDTpref(isfinite(bootDiffTDTpref)));
    TST.bootDiffTDTprefacc = 100*sum(isfinite(bootDiffTDTpref))/length(bootDiffTDTpref);

    % Anova based
    TST.bootTDT05    = median(bootTDT05(isfinite(bootTDT05)));
    TST.bootTDT05iqr =    iqr(bootTDT05(isfinite(bootTDT05)));
    TST.bootTDT05acc = 100*sum(isfinite(bootTDT05))/length(bootTDT05);

    TST.TDT01       = TDT01;
    TST.TDT01ci     = TDT01conf;

    TST.bootTDT01    = median(bootTDT01(isfinite(bootTDT01)));
    TST.bootTDT01iqr =    iqr(bootTDT01(isfinite(bootTDT01)));
    TST.bootTDT01acc = 100*sum(isfinite(bootTDT01))/length(bootTDT01);

    TST.TDT001      = TDT001;
    TST.TDT001ci    = TDT001conf;

    TST.bootTDT001    = median(bootTDT001(isfinite(bootTDT001)));
    TST.bootTDT001iqr =    iqr(bootTDT001(isfinite(bootTDT001)));
    TST.bootTDT001acc = 100*sum(isfinite(bootTDT001))/length(bootTDT001);

    % maximum response based
    TST.TDTmax05    = TDTmax05;
    TST.TDTmax05ci  = TDTmax05conf;

    TST.bootTDTmax05    = median(bootTDTmax05(isfinite(bootTDTmax05)));
    TST.bootTDTmax05iqr =    iqr(bootTDTmax05(isfinite(bootTDTmax05)));
    TST.bootTDTmax05acc = 100*sum(isfinite(bootTDTmax05))/length(bootTDTmax05);

    TST.TDTmax01    = TDTmax01;
    TST.TDTmax01ci  = TDTmax01conf;

    TST.bootTDTmax01    = median(bootTDTmax01(isfinite(bootTDTmax01)));
    TST.bootTDTmax01iqr =    iqr(bootTDTmax01(isfinite(bootTDTmax01)));
    TST.bootTDTmax01acc = 100*sum(isfinite(bootTDTmax01))/length(bootTDTmax01);

    TST.TDTmax001   = TDTmax001;
    TST.TDTmax001ci = TDTmax001conf;

    TST.bootTDTmax001    = median(bootTDTmax001(isfinite(bootTDTmax001)));
    TST.bootTDTmax001iqr =    iqr(bootTDTmax001(isfinite(bootTDTmax001)));
    TST.bootTDTmax001acc = 100*sum(isfinite(bootTDTmax001))/length(bootTDTmax001);

    % RF based
    TST.TDTrf05    = TDTrf05;
    TST.TDTrf05ci  = TDTrf05conf;

    TST.bootTDTrf05    = median(bootTDTrf05(isfinite(bootTDTrf05)));
    TST.bootTDTrf05iqr =    iqr(bootTDTrf05(isfinite(bootTDTrf05)));
    TST.bootTDTrf05acc = 100*sum(isfinite(bootTDTrf05))/length(bootTDTrf05);

    TST.TDTrf01    = TDTrf01;
    TST.TDTrf01ci  = TDTrf01conf;

    TST.bootTDTrf01    = median(bootTDTrf01(isfinite(bootTDTrf01)));
    TST.bootTDTrf01iqr =    iqr(bootTDTrf01(isfinite(bootTDTrf01)));
    TST.bootTDTrf01acc = 100*sum(isfinite(bootTDTrf01))/length(bootTDTrf01);

    TST.TDTrf001   = TDTrf001;
    TST.TDTrf001ci = TDTrf001conf;

    TST.bootTDTrf001    = median(bootTDTrf001(isfinite(bootTDTrf001)));
    TST.bootTDTrf001iqr =    iqr(bootTDTrf001(isfinite(bootTDTrf001)));
    TST.bootTDTrf001acc = 100*sum(isfinite(bootTDTrf001))/length(bootTDTrf001);


%     TST.DiffTDTpref     = DiffTDTpref;
%     TST.bootDiffTDTpref = bootDiffTDTpref;
%
%     TST.DiffTDTrf       = DiffTDTrf;
%     TST.bootDiffTDTrf   = bootDiffTDTrf;
%
%     TST.mnAng       = mnAng;
%     TST.mnAngDev    = mnAngDev;
%     TST.mnAngStd    = mnAngStd;
%
%     TST.mean_tune   = mean_tune;
%     TST.std_tune    = std_tune;
%     TST.AngMean_vec = AngMean;
%     TST.AngDev_vec  = AngDev;
%     TST.AngStd_vec  = AngStd;
%     TST.medianP     = median_pvec;
%     TST.upperP      = upperP;
%     TST.lowerP      = lowerP;
%
%     TST.medianPmax  = median_pMax;
%     TST.upperPmax   = upperPMax;
%     TST.lowerPmax   = lowerPMax;
end

%  ========================================================================
%% plot overview
if(doplot == 1)
    figure('Position', [0 0 1200 1000], 'Renderer', 'Painters');

    %% define plot colors
    col_sp1_main  = [  9  43   0]./255;
    col_sp1_area  = [134 236 107]./255;
    col_sp1_2nd   = [ 11  55   0]./255;

    col_sp2_main  = [  1   5  45]./255;
    col_sp2_area  = [109 121 235]./255;
    col_sp2_2nd   = [  1   5  45]./255;

    col_diff_main = [ 10  10  10]./255;
    col_diff_area = [128 128 128]./255;
    col_diff_2nd  = [ 25  25  25]./255;

    col_sig_main  = [255 125   0]./255;
    col_sig_area  = [255 190  56]./255;
    col_sig_2nd   = [255 125   0]./255;

    plotwin  = [ 0, 750];
    prefcnd = find(cndvec == mod(prefloc,8));
    anticnd = find(cndvec == mod(antiloc,8));

    %% response estimates without clipping
    subplot('Position',[0.05 0.76  0.40 0.2]);
    hold on;
    ylabel('activity [spikes/s]','fontsize',14)
    set(gca,'TickDir','out');
    set(gca,'fontsize',14);
    set(gca,'XTickLabel',[]);
    title(nfo, 'FontSize', 16);

    ppos = respraw.xtime >= plotwin(1) & respraw.xtime <= plotwin(2);

    
    for(i=1:length(cndvec))
        if(i == prefcnd || i == anticnd)
            continue;
        end
        vpos = isfinite(meanraw(i,:)) & isfinite(errraw(i,:)) & ppos == 1;
        patch([respraw.xtime(vpos), fliplr(respraw.xtime(vpos))],[meanraw(i,vpos)-errraw(i,vpos), fliplr(meanraw(i,vpos)+errraw(i,vpos))], [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
        plot(respraw.xtime(vpos), meanraw(i,vpos), 'k', 'linewidth', 2.5);
    end

    vpos = isfinite(meanraw(anticnd,:)) & isfinite(errraw(anticnd,:)) & ppos == 1;
    patch([respraw.xtime(vpos), fliplr(respraw.xtime(vpos))],[meanraw(anticnd,vpos)-errraw(anticnd,vpos), fliplr(meanraw(anticnd,vpos)+errraw(anticnd,vpos))], ...
        col_sp2_area, 'EdgeColor', col_sp2_2nd);
    plot(respraw.xtime(vpos), meanraw(anticnd,vpos), 'Color', col_sp2_main, 'linewidth', 2.5);

    vpos = isfinite(meanraw(prefcnd,:)) & isfinite(errraw(prefcnd,:)) & ppos == 1;
    patch([respraw.xtime(vpos), fliplr(respraw.xtime(vpos))],[meanraw(prefcnd,vpos)-errraw(prefcnd,vpos), fliplr(meanraw(prefcnd,vpos)+errraw(prefcnd,vpos))], ...
        col_sp1_area, 'EdgeColor', col_sp1_2nd);
    plot(respraw.xtime(vpos), meanraw(prefcnd,vpos), 'Color', col_sp1_main, 'linewidth', 2.5);
    vline(timwin,'color','g', 'linestyle', '--', 'linewidth',2.5);

    xlim(plotwin);

    %% raster overview
    subplot('Position',[0.05 0.05  0.40 0.7]);
    hold on;
    % title(fignm)
    xlabel('time [ms]','fontsize',14)
    %ylabel('condition blocks')
    set(gca,'TickDir','out');
    set(gca,'fontsize',14);
    set(gca, 'YTick',[], 'YTickLabel','')

    cnt = 0;
    for(i=0:length(cndvec)-1)
        currLoc = mod(prefloc+3+i,length(cndvec))+1;
        cpos = find(condvec == currLoc);

        cnd_spks = spks(cpos,:);
        cnd_spks(cnd_spks < plotwin(1) | cnd_spks > plotwin(2)) = NaN;

        if(isempty(cnd_spks))
            continue
        end
        
        cclip = cliptime(cpos);
        [~, cliporder] = sort(cclip);

        t1 = cnt;
        for(s=1:cndcnt(currLoc))
            cnt  = cnt + 1;

            tpos = cliporder(s);
            cspikes = cnd_spks(tpos,:);
            cspikes(isnan(cspikes)) = [];
            if(~isempty(cspikes))
                if(length(cliporder) == cndcnt(currLoc) && ~isempty(cliporder))
                    plot([cclip(tpos), cclip(tpos)],[cnt-0.5, cnt+0.5],'-r', 'LineWidth',1.5);
                end
                plot(cspikes,cnt,'.', 'color', col_sp1_2nd);
            end
        end

        plot([median(cclip), median(cclip)],[t1,cnt],'-r', 'LineWidth',1.5);

        cnt = cnt + 1;
        if(i == 3 || i == 4)
            hline(cnt,'Color',col_sp1_area, 'LineWidth',2);
            if(i == 4)
                plot(plotwin([1,1])+1,[t1,cnt],'color',col_sp1_area, 'LineWidth',4);
                plot(plotwin([2,2]),[t1,cnt],'color',col_sp1_area, 'LineWidth',4);
            end
        elseif(i==0)
            hline(cnt,'Color',col_sp2_area, 'LineWidth',2);
            hline(0,'Color',col_sp2_area, 'LineWidth',2);
            plot(plotwin([1,1])+1,[t1,cnt],'color',col_sp2_area, 'LineWidth',4);
            plot(plotwin([2,2]),[t1,cnt],'color',col_sp2_area, 'LineWidth',4);
        elseif(i<7)
            hline(cnt,'Color','red', 'LineWidth',1);
        end

        t = text(plotwin(1)-0.075*diff(plotwin),t1 + cndcnt(currLoc)/2, {['Loc ',int2str(currLoc-1)], [int2str(theta_vals(currLoc)),' deg']}, 'FontSize', 12);
        set(t,'HorizontalAlignment','center','VerticalAlignment','top', 'Rotation',90);
    end
    vline(timwin,'color','g', 'linestyle', '--', 'linewidth',2.5);

    ylim([0,cnt]);
    xlim(plotwin);

    %% response estimates after clipping
    subplot('Position',[0.5 0.68 0.48 0.275]);
    hold on;
    ylabel('activity [spikes/s]','fontsize',14)
    set(gca,'TickDir','out','fontsize',14);
    set(gca,'fontsize',14);
    set(gca,'XTickLabel',[]);

    for(i=1:length(cndvec))
        if(i == prefcnd || i == anticnd)
            continue;
        end
        vpos = isfinite(meanresp(i,:)) & isfinite(errresp(i,:)) & xpos == 1;

        patch([respest.xtime(vpos), fliplr(respest.xtime(vpos))],[meanresp(i,vpos)-errresp(i,vpos), fliplr(meanresp(i,vpos)+errresp(i,vpos))], [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
        plot(respest.xtime(vpos), meanresp(i,vpos), 'k', 'linewidth', 2.5);
    end

    vpos = isfinite(meanresp(anticnd,:)) & isfinite(errresp(anticnd,:)) & xpos == 1;
    patch([respest.xtime(vpos), fliplr(respest.xtime(vpos))],[meanresp(anticnd,vpos)-errresp(anticnd,vpos), fliplr(meanresp(anticnd,vpos)+errresp(anticnd,vpos))], ...
        col_sp2_area, 'EdgeColor', col_sp2_2nd);
    plot(respest.xtime(vpos), meanresp(anticnd,vpos), 'Color', col_sp2_main, 'linewidth', 2.5);

    vpos = isfinite(meanresp(prefcnd,:)) & isfinite(errresp(prefcnd,:)) & xpos == 1;
    patch([respest.xtime(vpos), fliplr(respest.xtime(vpos))],[meanresp(prefcnd,vpos)-errresp(prefcnd,vpos), fliplr(meanresp(prefcnd,vpos)+errresp(prefcnd,vpos))], ...
        col_sp1_area, 'EdgeColor', col_sp1_2nd);
    plot(respest.xtime(vpos), meanresp(prefcnd,vpos), 'Color', col_sp1_main, 'linewidth', 2.5);
    vline(timwin,'color','g', 'linestyle', '--', 'linewidth',2.5);

    vline(TDT05,    'color','r', 'linestyle',  '-', 'linewidth',2.5);
    vline(TDT05conf,'color','r', 'linestyle', '--', 'linewidth',2.5);

    xlim(xrng);

    %% response difference RF
    subplot('Position',[0.5 0.47 0.48 0.2]);
    hold on;
    ylabel('RF response difference','fontsize',14)
    set(gca,'TickDir','out');
    set(gca,'fontsize',14);
    set(gca,'XTickLabel',[]);

    vpos = isfinite(MaxBdiff) & respest.xtime > xrng(1) & respest.xtime < xrng(2);  %xpos == 1;
    if(~isempty(InRFpos) && ~isempty(OutRFpos) && length(MGRF) <= 4 && all(isfinite(MGRF)))
    %if(~isempty(MGRF) && all(isfinite(MGRF)))
        if(nBoot > 1)
            minval = (MaxBdiff(vpos) - MaxBstd(vpos));
            maxval = (MaxBdiff(vpos) + MaxBstd(vpos));
            patch([respest.xtime(vpos), fliplr(respest.xtime(vpos))], [minval, fliplr(maxval)], col_diff_area, 'EdgeColor', col_diff_2nd);
        end

        plot(respest.xtime(vpos), MaxBdiff(vpos), 'Color', col_diff_main, 'linewidth', 2.5);
        vline(timwin,'color','g', 'linestyle', '--', 'linewidth',2.5);

        vline(TST.bootDiffTDTrf,            'color','r', 'linestyle',  '-', 'linewidth',2.5);
        vline(prctile(bootDiffTDTrf,[5,95]),'color','r', 'linestyle', '--', 'linewidth',2.5);

        hline(prctile(MaxBdiff(spontpos),[1 99]), 'color','g', 'linestyle', '--', 'linewidth',2.5);
        hline(median(MaxBdiff(spontpos)),         'color','g', 'linestyle',  '-', 'linewidth',2.5);

        xlim(xrng);
   end

     %% response difference max
    subplot('Position',[0.5 0.26 0.48 0.2]);
    hold on;
    ylabel('Pref response difference','fontsize',14)
    set(gca,'TickDir','out');
    set(gca,'fontsize',14);
    set(gca,'XTickLabel',[]);

   if(~isempty(prefPos) && ~isempty(antiPos))
        if(nBoot > 1)
            minval = (PrefBdiff - PrefBstd);
            maxval = (PrefBdiff + PrefBstd);
            vpos = isfinite(minval) & isfinite(maxval);
            patch([respest.xtime(vpos), fliplr(respest.xtime(vpos))],[minval(vpos), fliplr(maxval(vpos))], col_diff_area, 'EdgeColor', col_diff_2nd);
        end

        plot(respest.xtime, PrefBdiff, 'Color', col_diff_main, 'linewidth', 2.5);
        vline(timwin,'color','g', 'linestyle', '--',  'linewidth',2.5);

        vline(TST.bootDiffTDTpref,            'color','r', 'linestyle',  '-', 'linewidth',2.5);
        vline(prctile(bootDiffTDTpref,[5,95]),'color','r', 'linestyle', '--', 'linewidth',2.5);

        hline(prctile(PrefBdiff(spontpos),[1 99]), 'color','g', 'linestyle', '--', 'linewidth',2.5);
        hline(median(PrefBdiff(spontpos)),         'color','g', 'linestyle',  '-', 'linewidth',2.5);

        xlim(xrng);
   end

    %% pvalue
    subplot('Position',[0.5 0.05 0.48 0.2]);
    hold on;
    ylabel('-log(p)','fontsize',14)
    xlabel('time [ms]','fontsize',14)
    set(gca,'TickDir','out');
    set(gca,'fontsize',14);

    if(nBoot > 1)
        %     minval = log(mean_pvec - std_pvec); maxval = log(mean_pvec +
        %     std_pvec);

        px = [xvec, fliplr(xvec)];
        py = -log([lowerP,    fliplr(upperP)]);
        px(isnan(py)) = [];
        py(isnan(py)) = [];
        patch(px, py,    col_sig_area,  'EdgeColor', col_sig_2nd);

        px = [xvec, fliplr(xvec)];
        py = -log([lowerPMax, fliplr(upperPMax)]);
        px(isnan(py)) = [];
        py(isnan(py)) = [];
        patch(px, py, col_diff_area, 'EdgeColor', col_diff_2nd);
    end

    plot(xvec, -log(median_pvec), 'Color', col_sig_main, 'linewidth', 2.5);
    plot(xvec, -log(median_pMax), 'Color', col_diff_main, 'linewidth', 2.5);

    p05  = hline(-log(0.05), 'color','red','LineStyle',':');
    p01  = hline(-log(0.01), 'color','red','LineStyle','--');
    p001 = hline(-log(0.001),'color','red','LineStyle','-');
    legend([p05 p01 p001], 'p = 0.05', 'p = 0.01', 'p = 0.001', 'Location', 'NorthEast');
    vline(timwin,'color','g', 'linestyle', '--', 'linewidth',2.5);
    %
    vline(TDT05,    'color','r', 'linestyle',  '-', 'linewidth',2.5);
    vline(TDT05conf,'color','r', 'linestyle', '--', 'linewidth',2.5);

    xlim(xrng);
    %ylim([-log(min(lowerP)), 0]);

    if(-log(min(lowerP)) > 0)
        ylim([0, -log(min(lowerP))]);
    end


%     %% tuning strength
%     subplot('Position',[0.5 0.47 0.48 0.2]);
%     hold on;
%     ylabel('tuning vector','fontsize',14)
%     set(gca,'TickDir','out');
%     set(gca,'fontsize',14);
%     set(gca,'XTickLabel',[]);
%
%     if(nBoot > 1)
%         minval = (mean_tune - std_tune);
%         maxval = (mean_tune + std_tune);
%         patch([xvec, fliplr(xvec)],[minval, fliplr(maxval)], col_diff_area, 'EdgeColor', col_diff_2nd);
%     end
%     plot(xvec, mean_tune, 'Color', col_diff_main, 'linewidth', 2.5);
%     vline(timwin,'color','g', 'linestyle', '--', 'linewidth',2.5);
%
%     vline(TDT05,    'color','r', 'linestyle',  '-', 'linewidth',2.5);
%     vline(TDT05conf,'color','r', 'linestyle', '--', 'linewidth',2.5);
%
%     xlim(xrng);
%
%     %% pref location
%     subplot('Position',[0.5 0.26 0.48 0.2]);
%     hold on;
%     ylabel('preferred location [deg]','fontsize',14)
%     set(gca,'TickDir','out');
%     set(gca,'fontsize',14);
%     set(gca,'XTickLabel',[]);
%
%     if(nBoot > 1)
%         minval = (AngMean - AngDev);
%         maxval = (AngMean + AngDev);
%         patch([xvec, fliplr(xvec)],[minval, fliplr(maxval)], col_diff_area, 'EdgeColor', col_diff_2nd);
%     end
%     plot(xvec, AngMean, 'Color', col_diff_main, 'linewidth', 2.5);
%
%     ylim([mnAng-180, mnAng+180]);
%     hline(mnAng,'color','r')
%     vline(timwin,'color','g', 'linestyle', '--', 'linewidth',2.5);
%
%     vline(TDT05,    'color','r', 'linestyle',  '-', 'linewidth',2.5);
%     vline(TDT05conf,'color','r', 'linestyle', '--', 'linewidth',2.5);
%
%     xlim(xrng);

end



%  ========================================================================
%% add nested functions

function TDTvec = TDT_getsigP(dtvec, xdata, thrval, minseq)

    sigmat = dtvec <= thrval;

    TDTvec = nan(size(dtvec,1),1);

    parfor(r=1:size(dtvec,1))
        sigpos   = min(strfind(sigmat(r,:), ones(1,minseq)));
        if(~isempty(sigpos))
            TDTvec(r) = xdata(sigpos);
        end
    end


%  ========================================================================

function [TDTest, TDTci] = TDT_getsigdiff(dtvec, uBnd, lBnd, xdata, thrval, minseq)

    TDTci = [NaN, NaN];

    if(~exist('xdata','var') || isempty(xdata))
        xdata = 1:length(dtvec);
    end

    issig  = dtvec <= thrval;

    TDTest = xdata(min(strfind(issig, ones(1,minseq))));
    if(isempty(TDTest))
        TDTest = NaN;
    end

    if(exist('lBnd','var') && ~isempty(lBnd))
        if(isfinite(TDTest))
            p = xdata(find(lBnd > thrval & xdata < TDTest,1,'last')+1);
            if(~isempty(p))
                TDTci(1) = p;
            end
        end
    end

    if(exist('uBnd','var') && ~isempty(uBnd))
        if(isfinite(TDTest))
            p = xdata(find(uBnd < thrval & xdata > TDTest,1,'first'));
            if(~isempty(p))
                TDTci(2) = p;
            end
        end
   end

