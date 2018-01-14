function SPKobj = wz_get_SPKobj(spk_times, tmwin, tzero, clip, na_val, prec, reso)
% SPKobj - Create a spike object (SPKobj).
%  
% DESCRIPTION 
% This function reates a structure summarizing spiking activity based on
% spike times. It's purpose is to simplify handling in other routines by
% keeping several relevant descriptors of the data in a single variable.
% 
% SYNTAX 
%   SPKobj = wz_get_SPKobj(spk_times, tmwin, tzero, clip, na_val, prec, reso)
%
% INPUT
%         <spk_times> should be a 2D matrix containing spike timies with rows representing
%                     trials. If spike times of 0 occurr they will be ignored (replaced by NaN).
% 
%         <tmwin>     determines the time window for spikes to be considered.
%                     default: [min(spk_times(:)) , max(spk_times(:))]
% 
%         <tzero>     specifies the time that should be set to zero. This could be a single scalar value or a
%                     vector with a length corresponding to the numbers of row (i.e. trials) of <spk_times>.
% 
%         <clip>      discard all spike times beyond this time. This could be a single scalar value or a
%                     vector with a length corresponding to the numbers of row (i.e. trials) of <spk_times>.
%                     Clipping is done before correcting times with <tzero>.
% 
%         <na_val>    is the identifier used for invalid spike entries that should be replaced by NaN.
% 
%         <prec>      sets the precision of the spike time. As default it is assumed that spike times
%                     are specified as 'ms'. In this case prec would be 1 (default).  If times are
%                     defined as 'sec' the precision would be 1000. Output times are specified in 'ms'.
% 
%         <reso>      temporal resolution of spike trains (and sampling of spike times). 
%                     Default is 1 [ms]. <reso> specifies the bin width of the spike train in ms,
%                     thus a resolution of 0.5 micro-seconds would be achieved with reso=0.5.
%                     Though <prec> and <reso> relate to each other, they are independent specifications. 
%                     It is, for example, possible to downsample data to a 1 ms resolution if the 
%                     input has a  precision of 0.1.
%            
% REFERENCES 
%
% 
%
% ......................................................................... 
% wolf zinke, wolfzinke@gmail.com 
%
% $Created : 11-Jan-2014 by wolf zinke
% $Modified: 31-Jul-2014 by wolf zinke

%% derive information from input
if(isfield(spk_times,'spk_times'))
    spk_times = spk_times.spk_times;
elseif(isfield(spk_times,'spiketimes'))
    spk_times = spk_times.spiketimes;  
end

SPKobj.nTrials  = size(spk_times,1);

% ========================================================================= 
%% define default values
% if(~exist('binwdth','var') || isempty(binwdth))
%    binwdth = 10;
% end

if(~exist('prec','var') || isempty(prec))
   prec = 1;
end

if(~exist('clip','var'))
   clip = [];
end

if(~exist('reso','var') || isempty(reso))
   reso = 1;
end

if(~exist('na_val','var'))
   na_val = []; % replace invalid entries with NaN
end

% ========================================================================= 
%% convert spike times to ms
spk_times = spk_times .* prec;

if(size(spk_times,1) > 1)
    spk_times = spk_times';
end

% ========================================================================= 
%% correct spike times by aligning them to an event specified by tzero.
if(exist('tzero','var') && ~isempty(tzero))
    if(~isempty(spk_times))
        spk_times = bsxfun(@minus, spk_times, tzero(:));
    else
        spk_times = NaN * tzero(:);
    end
    
    if(~isempty(clip))
        SPKobj.clip = bsxfun(@minus, clip(:), tzero(:)); % correct the clipping times to be aligned as well
    end
end

SPKobj.nTrials  = size(spk_times,1);
SPKobj.trialcnt = nan(SPKobj.nTrials,1);

% ========================================================================= 
%% clip spike times
if(exist('clip','var') && ~isempty(clip))
    if(length(clip) == 1)
        % spk_times(spk_times>clip) =  NaN;
        clip = repmat(clip,SPKobj.nTrials,1);
        SPKobj.clip = clip;
        SPKobj.trial_order = 1:SPKobj.nTrials;
    elseif(~isempty(spk_times))
        % spk_times(bsxfun(@gt, spk_times, clip(:))) = NaN;
        [~ , ordr] = sort(clip,'descend');
        SPKobj.clip = clip;
        SPKobj.trial_order = ordr;
    end
else
    SPKobj.clip = [];
    SPKobj.trial_order = 1:SPKobj.nTrials;
end

SPKobj.clipped = false;  % just keep information for clipping, but do not do that here.

% ========================================================================= 
%% if not specified derive relevant time window from data
if(~exist('tmwin','var') || isempty(tmwin))
   tmwin = [floor(min(spk_times(:))) , ceil(max(spk_times(:)))];
end

% ========================================================================= 
%% define values that should be considered as NaN.
if(exist('na_val','var') && ~isempty(na_val))
   spk_times(spk_times == na_val) = NaN; % replace invalid entries with NaN
end

% ========================================================================= 
%% constrain valid spikes to a given time window and crop NaN values
SPKobj.spiketimes = wz_cropNaN(spk_times, tmwin(1), tmwin(2), na_val);
SPKobj.timewindow = tmwin;

%% get number of spikes for each trial
for(t=1:SPKobj.nTrials)
    SPKobj.trialcnt(t) = sum(~isnan(spk_times(t,:)));
end

% ========================================================================= 
%% get the spike time matrix and spike train matrix
SPKobj.resolution = reso;
SPKobj.traintime  = tmwin(1)-reso : reso : tmwin(2);
SPKobj.spiketrain = nan(SPKobj.nTrials, length(SPKobj.traintime));

for(t=1:SPKobj.nTrials)
    cSpk = spk_times(t,:);
    cSpk(isnan(cSpk)) = [];
    cSpk(cSpk < SPKobj.traintime(1) | cSpk > SPKobj.traintime(end)) = [];
    
    SPKobj.spiketrain(t,:) = hist(cSpk, SPKobj.traintime);
end

% ========================================================================= 
%% generate a spike histogram
% if(~isempty(tmwin))
%     % get spike histogram
%     xhist = tmwin(1)-binwdth/2 : binwdth : tmwin(2)+binwdth/2;
%     SPKobj.hist = [xhist + binwdth/2 ; histc(SPKobj.spiketimes(:), xhist)'];
%     SPKobj.hist(2,:) = (SPKobj.hist(2,:) ./ SPKobj.nTrials) * (1000/binwdth); % first term gets the trial average, second one converts spikes/s.
%     % bar(SPKobj.hist(1,:), SPKobj.hist(2,:));
%     SPKobj.histbin = binwdth;
% 
%     % get the cumulative sum
%     spkhist = histc(spk_times(:), [tmwin(1)-0.5 : tmwin(2)]+0.5);
%     SPKobj.cumsum = cumsum(spkhist)';
%     SPKobj.cumsum = SPKobj.cumsum ./ SPKobj.cumsum(end);
% else
%     SPKobj.hist   = [];
%     SPKobj.cumsum = [];
% end
% 
% 

