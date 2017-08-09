function Tdiff = wz_spk_GetRespDiff(, timwin, mincons, pthr, nBoot, smplsz)
% wz_spk_GetRespDiff - this is a wrapper function that prepared data for example
%                      for wz_spk_RespDiff_ranksum and aims for the identification
%                      of the time where two response profiles clearly diverge 
%                      
%  
% DESCRIPTION 
% 	
%  
% SYNTAX 
%   TST = wz_spk_GetRespDiff()
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
% $Created : 08-Oct-2014 by wolf zinke
% $Modified: 

%  ========================================================================
%% check input arguments and define defaults

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

if(~exist('timwin','var') || isempty(timwin))
    timwin = [1, min(size(cnd1, 1), size(cnd2, 1))];
end

%  ========================================================================
%% Define constant parameter

theta_vals = [0:45:315];
% thv        = theta_vals ./ (180/pi);
% angstep    = 45/(180/pi);


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

%  ========================================================================
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
    smplsz = min(cndcnt);
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

prefPos = find(condvec == prefloc-1);
antiPos = find(condvec == antiloc-1);

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

