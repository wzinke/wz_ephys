function Tdiff = wz_spk_RespDiff_ranksum(cnd1, cnd2, timwin, mincons, pthr, nBoot, smplsz, xvec)
% wz_spk_RespDiff_ranksum - calculate the time when two responses become
% clearly ifferen
%  
% DESCRIPTION 
% 	inspired by the function getTDT_SP provided by Rich Heitz
%  
% SYNTAX 
%   TST = wz_spk_RespDiff_ranksum()
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

% make sure both matrices have rows of the same length
if(size(cnd1,2) ~= size(cnd2,2))
    if(size(cnd1,2) > size(cnd2,2))
        cnd1(:,size(cnd2,2):end) = [];
    else
        cnd2(:,size(cnd1,2):end) = [];
    end
end

if(~exist('xvec','var') || isempty(xvec))
    xvec = 1:size(cnd1,2);
end

if(~exist('timwin','var') || isempty(timwin))
    timwin = [min(xvec), max(xvec)];
end

% =========================================================================
%% check and prepare data
Ntrials_1 = size(cnd1, 1);
Ntrials_2 = size(cnd2, 1);

Nbins = diff(timwin)+1;

% if there is no data just do nothing
if(Ntrials_1 < 4 || Ntrials_2 < 4)
    warning('Not sufficient trials -> aborted')
    Tdiff = [];
    return;
end

% matrix for bootstrap based bin-wise p-values
pMat = nan(nBoot, Nbins);

% =========================================================================
%% create matrix for bootstrap repetitions    
 
if(nBoot > 1)
    rng(sum(100*clock), 'twister');
    
    if(isempty(smplsz))
        smplsz1 = Ntrials_1;
        smplsz2 = Ntrials_2;
    else
        smplsz1 = smplsz;
        smplsz2 = smplsz;
    end
    
    Bsmpl1 = randi(Ntrials_1, nBoot, smplsz1);
    Bsmpl2 = randi(Ntrials_2, nBoot, smplsz2);
end

% make sure that the first sample is the 'raw data'
if(isempty(smplsz) || nBoot == 1)
    Bsmpl1(1,:) = 1:Ntrials_1;
    Bsmpl2(1,:) = 1:Ntrials_2;
end

% =========================================================================
%% go through each time bin and run a rank sum test on trial based firing rates
tpos = find(xvec >= timwin(1) & xvec <= timwin(2));

parfor(i=1:length(tpos))  % loop over time bins
    cBinVal1 = cnd1(:, tpos(i));
    cBinVal2 = cnd2(:, tpos(i));

    for(b=1:nBoot) % loop over bootstrap repetitions
        cBSmpl1 = cBinVal1(Bsmpl1(b,:)); % resample data
        cBSmpl1(isnan(cBSmpl1)) = [];    % eliminate nan values
        cBSmpl2 = cBinVal2(Bsmpl2(b,:));
        cBSmpl2(isnan(cBSmpl2)) = [];

        if(length(cBSmpl1) >= 4 && length(cBSmpl2) >= 4)
            pMat(b,i) = ranksum(cBSmpl1, cBSmpl2, 'method','approximate');
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
Tdiff.Resp1 = nanmean(cnd1);
tcnt = sum(isfinite(cnd1));
Tdiff.Std1 = nanstd(cnd1) ./ sqrt(tcnt);

Tdiff.Resp2 = nanmean(cnd2);
tcnt = sum(isfinite(cnd2));
Tdiff.Std2 = nanstd(cnd2) ./ sqrt(tcnt);

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




