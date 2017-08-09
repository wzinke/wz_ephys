function Tdiff = wz_spk_RespTune(respmat, condvec, timwin, mincons, pthr, nBoot, smplsz, xvec)
% wz_spk_RespTune - determine tuning dynamics
%
%  
% DESCRIPTION 
%
%  
% SYNTAX 
%   TST = wz_spk_RespTune()
%
%   Input:
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

%  ========================================================================
%% check input arguments and define defaults
if(nargin < 2)
    error('At least two arguments required!');
end

% number consecutive bins exceeding the threshold
if(~exist('mincons','var') || isempty(mincons))
    mincons = 20;
end

% p-value threshold
if(~exist('pthr','var') || isempty(pthr))
    pthr = 0.01;
end

% number of bootstrap repetitions
if(~exist('nBoot','var') || isempty(nBoot))
    nBoot = 1;
end

% number of bootstrap repetitions
if(~exist('smplsz','var') || isempty(smplsz))
    smplsz = [];
end

if(~exist('xvec','var') || isempty(xvec))
    xvec = 1:size(respmat,2);
end

if(~exist('timwin','var') || isempty(timwin))
    timwin = [min(xvec), max(xvec)];
end

% =========================================================================
%% check and prepare data
Nbins = diff(timwin)+1;

% if there is no data just do nothing
if(Ntrials_1 < 4 || Ntrials_2 < 4)
    warning('Not sufficient trials -> aborted')
    Tdiff = [];
    return;
end

% matrix for bootstrap based bin-wise p-values
pMat = nan(nBoot, Nbins);

% get list of available conditions
[cndcnt, cndvec] = get_elecounts(condvec);
% cndvec = unique(condvec);
Ncond   = length(cndvec);
% cndcnt = hist(condvec, 1:Ncond); % how many trials per condition?

if(any(cndcnt < 4))
    warning('Conditions with less than 4 trials deteted. They will be removed from analysis!')
    
    lesscond = cndvec(cndcnt < 4);
    rempos = ismember(condvec,lesscond);
    
    respmat(rempos,:)  = [];
    condvec(rempos)    = [];
%     cndvec(cndcnt < 4) = [];
%     Ncond  = length(cndvec);
%     cndcnt = hist(condvec, 1:Ncond);
    [cndcnt, cndvec] = get_elecounts(condvec);
    Ncond  = length(cndvec);
end

Ntrials =    sum(cndcnt);

if(Ncond < 2)
    error('Only one condition - nothing to compare!')
end

% =========================================================================
%% create matrix for bootstrap repetitions
if(isempty(smplsz))
    Bsmpl   = nan(nBoot, Ntrials);
    smplszs = cndcnt;
else
    Bsmpl   = nan(nBoot,smplsz*length(cndvec));
    smplszs = rep(smplsz, Ncond,1);
end
    
if(nBoot > 1)
    rng(sum(100*clock), 'twister');
     
    for(i=1:Ncond)
        cpos = find(condvec == cndvec(i));
        if(i==1)
            cs = 1;
        else
            cs = cumsum(smplszs(i-1));
        end
        
        Bsmpl(:,cs:cs+smplszs(i)-1) = cpos(randi(length(cpos), nBoot, smplszs(i)));
    end
end

% make sure that the first sample is the 'raw data'
if(isempty(smplsz) || nBoot == 1)
    for(i=1:Ncond)
        cpos = find(condvec == cndvec(i));
        if(i==1)
            cs = 1;
        else
            cs = cumsum(smplszs(i-1));
        end        
        Bsmpl(1,cs:cs+smplszs(i)-1) = cpos;
    end
end

% =========================================================================
%% go through each time bin and run a rank sum test on trial based firing rates
tpos = find(xvec >= timwin(1) & xvec <= timwin(2));

parfor(i=1:length(tpos))  % loop over time bins
    cBinVals = respmat(:, tpos(i));

    for(b=1:nBoot) % loop over bootstrap repetitions
        csr  = cBinVals(Bsmpl(b,:),i); % resample data
        cloc = condvec( Bsmpl(b,:));
        cloc(isnan(csr)) = [];         % eliminate nan values
        csr( isnan(csr)) = [];

        ccndcnt = hist(cloc, 1:length(unique(cloc))); % how many trials per condition?

        if(any(ccndcnt < 4))
            lesscond = cndvec(cndcnt < 4);
            rempos   = ismember(cloc,lesscond);

            csr(rempos)  = [];
            cloc(rempos) = [];
        end
        
        if(length(unique(cloc)) >= 2)
            pMat(b,i) = kruskalwallis(csr, cloc, 'off');
        end
    end
end

% =========================================================================
%%  determine time point of clear response difference
if(nBoot > 1)
    median_Tdiff = nanmedian(pMat,1);
    upper_Tdiff  =   prctile(pMat, 75, 1);
    lower_Tdiff  =   prctile(pMat, 25, 1);

    mean_Tdiff = nanmean(pMat,1);
    std_Tdiff  = nanstd(pMat,1);
end

Tdiff_bootvec = nan(nBoot, length(pthr));
for(p=1:length(pthr))
    if(nBoot > 1)
        [TDTmed(p),   TDTmedCI(p,:)] = wz_spk_getTSdiff(median_Tdiff, pthr(p), mincons, upper_Tdiff, lower_Tdiff,xvec(tpos));
        [TDTmean(p), TDTmeanCI(p,:)] = wz_spk_getTSdiff(mean_Tdiff,   pthr(p), mincons, mean_Tdiff+std_Tdiff, mean_Tdiff-std_Tdiff,xvec(tpos));
    end
    
    Tdiff_bootvec(:,p) = wz_spk_getTSdiff(pMat, pthr(p), mincons,[],[],xvec(tpos));
end

% =========================================================================
%% create output struct
cndvec = unique(condvec);
Ncond  =  length(cndvec);

Tdiff.Ncond = Ncond;
Tdiff.trials_per_cond = sprintf('%d_',cndcnt);

for(i=1:Ncond)
    cpos = find(condvec == cndvec(i));
    Tdiff.Resp(i,:) = nanmean(respmat(cpos,:));
    tcnt = sum(isfinite(respmat(cpos,:)));
    Tdiff.Std1 = nanstd(respmat(cpos,:)) ./ sqrt(tcnt);
end

Tdiff.xvec = xvec;

Tdiff.pMat = pMat;

if(nBoot > 1)
    if(isempty(smplsz))
        Tdiff.DiffTime_raw = Tdiff_bootvec(1);
    end
    
    Tdiff.DiffTime     = nanmedian(Tdiff_bootvec);
    Tdiff.DiffTime_CI  = prctile(Tdiff_bootvec, [25 75], 1);
    Tdiff.DiffTime_acc = sum(isfinite(Tdiff_bootvec)) ./ nBoot;
    
    Tdiff.Tdiff_boot = Tdiff_bootvec;
    
    Tdiff.TDTmed    = TDTmed;
    Tdiff.TDTmedCI  = TDTmedCI;
    
    Tdiff.TDTmean   = TDTmean;
    Tdiff.TDTmeanCI = TDTmeanCI;
else
    Tdiff.DiffTime = Tdiff_bootvec;
end




