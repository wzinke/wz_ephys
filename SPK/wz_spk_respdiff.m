function [TDTest , TDT] = wz_spk_respdiff(spk1, spk2, clip1, clip2, timwin, nBoot, do_plot, fignm)
% Determines the time where two response profiles start to become different
% from each other. 
%
% wolf zinke, 14.1.2014
% TODO: limit number of trials used for calculations to be the same across
%       conditions, or limit it to a specified trial number.

siglev  = 0.01;   % Wilcoxon ranksum test must give a bin wise p-value below that level.
mincons = 10;     % at least <mincons> time points below the <siglev> to be considered

%% check input and define default values
if(~isfield(spk1,'spikedensities'))
    spk1 = wz_spk_density(spk1, 'exp', 20);
end

if(~isfield(spk2,'spikedensities'))
    spk2 = wz_spk_density(spk2, 'exp', 20);
end

% ensure a time window is used for the test that suits both series 
if(isempty(spk1.timewindow) || isempty(spk2.timewindow))
    warning('At least one of the spike matrices do not contain any data! \n ... nothing to compare here.');
    TDTest = NaN;
    TDT = [];
    return
end

if(~exist('timwin','var') || isempty(timwin))
    timwin(1) = max([spk1.timewindow(1); spk2.timewindow(1)]);
    timwin(2) = min([spk1.timewindow(2); spk2.timewindow(2)]);
elseif(any([spk1.timewindow(1);spk2.timewindow(1)] > timwin(1)) || ...
       any([spk1.timewindow(2);spk2.timewindow(2)] < timwin(2)) )
    timwin(1) = max([spk1.timewindow(1); spk2.timewindow(1)]);
    timwin(2) = min([spk1.timewindow(2); spk2.timewindow(2)]);
end

% number of bootstrap repetitions
if(~exist('nBoot','var') || isempty(nBoot))
    nBoot = 1;
end

% minimum requirement of consecutive significant time bins
if(~exist('do_plot','var') || isempty(do_plot))
    do_plot = 0;
end

% minimum requirement of consecutive significant time bins
if(~exist('fignm','var') )
    fignm = [];
else
    do_plot = 1;
end

% % minimum requirement of consecutive significant time bins
% if(~exist('mincons','var') || isempty(mincons))
%     mincons = 10;
% end
% 
% % minimum requirement of consecutive significant time bins
% if(~exist('siglev','var') || isempty(siglev))
%     siglev = 0.01;
% end

%% get the trial based spike densities
tp1 = spk1.xtime >= timwin(1) & spk1.xtime <= timwin(2);
tp2 = spk2.xtime >= timwin(1) & spk2.xtime <= timwin(2);

tvec = timwin(1):timwin(2);  % time bins used for the test

clipped = [0 0];
% clip spike trains if required (i.e. set to NaN)
if(exist('clip1','var') && ~isempty(clip1))
    spk1b = wz_spk_density(spk1, spk1.kerneltype, spk1.kernelwidth, spk1.nBoot, clip1);

%     tmmat = repmat(tvec,spk1.nTrials,1);
%     spkmat1(bsxfun(@gt, tmmat, clip1(:))) = NaN;
%     [~ , ordr] = sort(clip1,'descend');
%     spk1.trial_order = ordr;
    
    maxt1 = find(sum(~isnan(spk1b.spikedensities(:, tp1))) < 4,1,'first'); % max reasonable length
    if(isempty(maxt1))
        maxt1 = length(tvec);
    end
    
    if(prctile(spk1b.clip,75) < timwin(2))
        timwin(2) = prctile(spk1b.clip,75);
    end
    clipped(1) = 1;
else
    spk1b = spk1; 
    maxt1 = length(tvec);
end

% clip spike trains if required (i.e. set to NaN)
if(exist('clip2','var') && ~isempty(clip2))
    spk2b = wz_spk_density(spk2, spk2.kerneltype, spk2.kernelwidth, spk2.nBoot, clip2);
%     tmmat = repmat(tvec,spk2.nTrials,1);
%     spkmat2(bsxfun(@gt, tmmat, clip2(:))) = NaN;
%     [~ , ordr] = sort(clip2,'descend');
%     spk2.trial_order = ordr;
    
    maxt2 = find(sum(~isnan(spk2b.spikedensities(:, tp2))) < 4,1,'first'); % max reasonable length
    if(isempty(maxt2))
        maxt2 = length(tvec);
    end
%     if(median(spk2b.clip) < timwin(2))
%         timwin(2) = median(spk2b.clip);
    if(prctile(spk2b.clip,75) < timwin(2))
        timwin(2) = prctile(spk2b.clip,75);
    end
    clipped(2) = 1;
else
    spk2b = spk2; 
    maxt2 = length(tvec);
end

% cut down analysis window to a reasonable length, 
% i.e. minimum sample size of 4 for both groups.
if( maxt1 < diff(timwin)+1 || maxt2 < diff(timwin)+1 )
    maxt = min([maxt1 maxt2]);
    
    tvec(maxt+1:end)      = [];
    timwin(2) = max(tvec);
else
    tvec = timwin(1):timwin(2); % correct 
end

tp1 = spk1b.xtime >= timwin(1) & spk1b.xtime <= timwin(2);
tp2 = spk2b.xtime >= timwin(1) & spk2b.xtime <= timwin(2);

spkmat1   = spk1b.spikedensities(:, tp1);
spkmean1  = nanmean(spkmat1,1);
trialcnt1 = sum(~isnan(spkmat1));
spkste1   = nanstd(spkmat1,1) ./ sqrt(trialcnt1);

spkmat2   = spk2b.spikedensities(:, tp2);
spkmean2  = nanmean(spkmat2,1);
trialcnt2 = sum(~isnan(spkmat2));
spkste2   = nanstd(spkmat2,1) ./ sqrt(trialcnt2);

if(timwin(1)==timwin(2))
    TDTest = NaN;
    TDT = [];
    return
end

%% pre-allocate matrices for the results
pvals = nan(nBoot,length(tvec));

% get the boot sample matrix
rng(sum(100*clock),'twister');

Bsmpl1 = randi(spk1b.nTrials, spk1b.nTrials, nBoot);
Bsmpl1(:,1) = 1:spk1b.nTrials; % make sure that first sample is original data
Bsmpl2 = randi(spk2b.nTrials, spk2b.nTrials, nBoot);
Bsmpl2(:,1) = 1:spk2b.nTrials;

bootmeans1 = nan(nBoot,size(spkmat1,2));
bootmeans2 = nan(nBoot,size(spkmat2,2));

%% run bootstrap loops 
for(r=1:nBoot)
    % get the current trial sample (with replacement)   
    cspikes1 = spkmat1(Bsmpl1(:,r),:);
    cspikes2 = spkmat2(Bsmpl2(:,r),:);
    
    bootmeans1(r,:) = nanmean(cspikes1,1);
    bootmeans2(r,:) = nanmean(cspikes2,1);

    maxtest = min([find(trialcnt1 < 4,1,'first'), find(trialcnt2 < 4,1,'first')])-tvec(1); % max reasonable length
    if(isempty(maxtest))
        maxtest = length(tvec);
    end
    
    parfor(b=1:maxtest)
        spk_one = cspikes1(:,b);
        spk_one(isnan(spk_one)) = [];
        spk_two = cspikes2(:,b);
        spk_two(isnan(spk_two)) = [];
        
        if(length(spk_one) > 4 && length(spk_two) > 4)
            pvals(r,b) = ranksum(cspikes1(:,b), cspikes2(:,b));  
        end
    end
end

%% identify first occurrence of the response differences
issig = pvals <= siglev;

TDTvals = nan(1,nBoot);
TDTp05  = nan(1,nBoot);
TDTp001 = nan(1,nBoot);
pdist   = nan(1,nBoot);

% there must be a faster and more elegant way instead of looping again!
for(r=1:nBoot) 
    sigpos = min(strfind(issig(r,:), ones(1,mincons)));
    
    if(~isempty(sigpos))
        TDTvals(r) = tvec(sigpos);
        
        pos1 = find(pvals(r,1:sigpos) > 0.05,1,'last');
        if(~isempty(pos1))
            pos = find(pvals(r,pos1:sigpos) < 0.05,1,'first');
            
            if(~isempty(pos))
                TDTp05(r) = tvec(pos+pos1); 
                
                pdist(r) = TDTvals(r) - TDTp05(r);
            end
        end
        
        pos = find(pvals(r,sigpos:end) < 0.001,1,'first');
        if(~isempty(pos))
            TDTp001(r) = tvec(sigpos+pos-1);            
        end
    end
end

% use the quartiles of the bootstrapped p-values to get confidence intervals
mP = nanmedian(pvals,1);
critpass = strfind(mP<siglev, ones(1,mincons));
critfail = strfind(mP>siglev, ones(1,mincons));

if(~isempty(critpass))
    thrpos = critpass(1);
    pTDT   = tvec(thrpos); 
else
    pTDT = NaN;
end

if(nargout ==2)
    upperP = prctile(pvals,75,1);
    pTDTu = tvec(min(strfind(upperP<siglev, ones(1,mincons)))); 
    if(isempty(pTDTu))
        pTDTu = NaN;
    end

    lowerP = prctile(pvals,25,1); 
    pTDTl = tvec(min(strfind(lowerP<siglev, ones(1,mincons)))); 
    if(isempty(pTDTl))
        pTDTl = NaN;
    end

    if(~isempty(critpass))
        % try some preliminary simple distance measures
        logpvec = abs(log(mP));
        % area below p-threshold
        sigpos = mP < siglev & tvec > pTDT;
    
        sigarea     = sum(logpvec(sigpos));
        sigarea_eff = sigarea/sum(sigpos);
    
        lp = critpass(1);
        sigarea_mat = [];
        while(~isempty(lp))
    %         np = find(mP >= siglev & tvec > tvec(lp),1,'first');
            np = min(critfail(critfail>lp));
            if(~isempty(np))
                carea = sum(logpvec(lp:np));
                sigarea_mat = [sigarea_mat; tvec([lp, np-1]) carea carea/length(lp:np)];
            else
                carea = sum(logpvec(lp:end));
                sigarea_mat = [sigarea_mat; tvec([lp, length(tvec)]) carea carea/length(lp:np)];
                break;
            end

            lp = min(critpass(critpass>np));
        end
    else
        sigarea = NaN;
        sigarea_eff = NaN;
        sigarea_mat = nan(1,4);
    end
end
%% look into the response difference
diffmat = bootmeans1 - bootmeans2;

%% generate output
TDTest = median(pTDT);

if(nargout == 2)
    TDT.spk1        = spk1;
    TDT.spk2        = spk2;
    TDT.times       = tvec;
    TDT.spkmat1     = spkmat1;
    TDT.spkmean1    = spkmean1;
    TDT.trialcnt1   = trialcnt1;
    TDT.spkste1     = spkste1;
    TDT.spkmat2     = spkmat2;
    TDT.spkmean2    = spkmean2;
    TDT.trialcnt2   = trialcnt2;
    TDT.spkste2     = spkste2;
    TDT.clipped     = clipped;
    TDT.timeWin     = timwin;
    TDT.meandiff    = nanmean(diffmat,1);
    TDT.stddiff     = nanstd(diffmat,0,1);
    TDT.mediandiff  = nanmedian(diffmat,1);
    TDT.upperdiff   = prctile(diffmat,75,1);
    TDT.lowerdiff   = prctile(diffmat,25,1);
    TDT.Pmat        = pvals;
    TDT.mP          = mP;
    TDT.uP          = upperP;
    TDT.lP          = lowerP;
    TDT.TDT         = TDTvals;
    TDT.acc         = sum(~isnan(TDTvals))/length(TDTvals);
    TDT.conf        = prctile(TDT.TDT,[5 95]);
    TDT.pTDT        = [pTDTl pTDT pTDTu];
    TDT.medianTDT   = nanmedian(TDTvals);
    TDT.meanTDT     = nanmean(  TDTvals);
    TDT.TDTp05      = TDTp05;
    TDT.TDTp001     = TDTp001;
    TDT.pdist       = pdist;
    TDT.medTDTp05   = nanmedian(TDTp05);
    TDT.medTDTp001  = nanmedian(TDTp001);
    TDT.medpdist    = nanmedian(pdist);
    TDT.sigarea     = sigarea;
    TDT.sigarea_eff = sigarea_eff;
    TDT.sigarea_mat = sigarea_mat;
end

%% create plots
if(do_plot ~= 0)
   wz_spk_plot_TDR(TDT, timwin, 1, fignm);
end

