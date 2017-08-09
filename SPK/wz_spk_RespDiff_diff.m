function Tdiff = wz_spk_RespDiff_diff(cnd1, cnd2, timwin, spontwin, mincons, nBoot, smplsz, xvec)
% wz_spk_RespDiff_diff - calculate the time when two responses become clearly different
%  
% DESCRIPTION 
% 	inspired by the function getTDT_SP provided by Rich Heitz
%  
% SYNTAX 
%   TST = wz_spk_RespDiff_diff()
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

if(~exist('spontwin','var') || isempty(spontwin))
    spontwin = [-500, 0];
end

% =========================================================================
%% check and prepare data
Ntrials_1 = size(cnd1, 1);
Ntrials_2 = size(cnd2, 1);

% if there is no data just do nothing
if(Ntrials_1 < 4 || Ntrials_2 < 4)
    warning('Not sufficient trials -> aborted')
    Tdiff = [];
    return;
end

% =========================================================================
%% get estimate based on response differences
tpos     = xvec >=   timwin(1) & xvec <=   timwin(2);
spontpos = xvec >= spontwin(1) & xvec <= spontwin(2);

RespDiff    = NaN; 
RespDiff_CI = [NaN, NaN];
bootDiff    = nan(1,nBoot);

[Bdiff, Bstd, bootdiff] = SPK_boot_RespDiff(cnd1, cnd2, nBoot, smplsz);

Bdiff    = abs(Bdiff);                                              
% DiffXced = Bdiff > prctile(Bdiff(spontpos)+Bstd(spontpos), 99) & tpos == 1;
DiffXced = Bdiff > median(Bdiff(spontpos))+3*mad(Bdiff(spontpos)) & tpos == 1;

sigpos   = min(strfind(DiffXced, ones(1,mincons)));

if(~isempty(sigpos))
    RespDiff = xvec(sigpos); 
    
%     DiffXced = Bdiff-Bstd > prctile(Bdiff(spontpos)+Bstd(spontpos), 99) & tpos == 1 & xvec < RespDiff;
    DiffXced = Bdiff-Bstd > median(Bdiff(spontpos))+3*mad(Bdiff(spontpos)) & tpos == 1 & xvec < RespDiff;
    lCI = min(strfind(DiffXced, ones(1,mincons)));
    if(~isempty(lCI))
        RespDiff_CI(2) = lCI;
    end
    
%     DiffXced = Bdiff+Bstd > prctile(Bdiff(spontpos)+Bstd(spontpos), 99) & tpos == 1 & xvec > RespDiff;
    DiffXced = Bdiff+-Bstd > median(Bdiff(spontpos))+3*mad(Bdiff(spontpos)) & tpos == 1 & xvec > RespDiff;
    uCI = min(strfind(DiffXced, ones(1,mincons)));
    if(~isempty(uCI))
        RespDiff_CI(2) = uCI;
    end
end

bootdiff = abs(bootdiff);

DiffXced_mat = bsxfun(@gt, bootdiff, prctile(bootdiff(:,spontpos), 99,2)) .* repmat(tpos,nBoot,1);

if(nBoot > 1)
    parfor(r=1:nBoot)
        
        sigpos = min(strfind(DiffXced_mat(r,:), ones(1,mincons)));
        if(~isempty(sigpos))
            bootDiff(r) = xvec(sigpos);
        end
    end
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

if(nBoot > 1)
    if(isempty(smplsz))
        Tdiff.DiffTime_raw    = bootDiff(1);
        Tdiff.DiffTime_raw_CI = RespDiff_CI;
    end
    
    Tdiff.DiffTime     = nanmedian(bootDiff);
    Tdiff.DiffTime_CI  = prctile(bootDiff, [25 75], 1);
    Tdiff.DiffTime_acc = sum(isfinite(bootDiff)) ./ nBoot;
    
    Tdiff.Tdiff_boot = bootDiff;
else
    Tdiff.DiffTime = RespDiff;
end




