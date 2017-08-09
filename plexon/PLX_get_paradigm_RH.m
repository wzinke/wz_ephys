function plxMat = PLX_get_paradigm_RH(plxMat, EV, rig, NHP, keep_sess, skip_waves, keep_nonstim)
% Read in a plexon data file organize it according to the paradigm information.
%
% DESCRIPTION
%
%
% SYNTAX
%
%   [plxdt] = PLX_get_paradigm(plxfile, EV, PostTrialTime, MedTrialSpike, SampleRate)
%
%   Input:
%
%       plxfile       - file name of the plexon data file (it also accept the
%                       struct that is generated with PLX2Mat - this might be
%                       helpful for debugging since it does not be read in
%                       every time).
%
%       EV            - A structure that defines event codes or a tempo EVENTDEF file
%
%       prdgm         - string that identify the paradigm and rig set up in
%                       order to load a default configuration.
%
%       ADchan        - Assignment of meanings of the analog channels.
%
%       keep_sess     - keep analoge data and spike times as single vector for
%                       the complete session (memory demand increases considerably)
%
%       skip_waves    - do not keep single spike wave forms after the average
%                       is determined. This will reduce memory usage a lot.
%
%       keep_nonstim  - do not discard trials with no stimulus onset
%
% .........................................................................
% wolf zinke, wolfzinke@gmail.com
%
% $Created : 26-Jan-2015 by wolf zinke
%
% ToDo:   - LFP channel identification is a bit ad hoc and might need some more
%           a priori information to work robustly.
%
%         - Get corrected SRT values (hack a function to determine the movement
%           onset of a saccade by first finding the movement distance exceeding
%           a threshold, and then identify the time or rise)
%
%         - outsource the event code interpretation to allow flexibility of the code
%
%         - determine ITI/ISI

% ____________________________________________________________________________%
%% Set some default values
PreTrialTime  = -500; % start sampling trial based data at this time relative to trial onset. MAKE IT LONG! Trial start defined as stimulus onset!
PostTrialTime =  250; % continue to collect analog data for this time period
PreOnsetTime  = -750; % discard analog data before this time, pad with nans if this time is not matched
SampleRate    = 1000;

% ____________________________________________________________________________%
%% get information about the current file version
% if(exist('wz_get_git_version','file'))
%     [gitversion, ~, GitInfo]   = wz_get_git_version(mfilename('fullpath'));
%     plxMat.ParadigmFileVersion = gitversion;
%     plxMat.ParadigmFileDate    = GitInfo.date;
% end

% ____________________________________________________________________________%
%% check input data
% if not specified use GUI to get the file
if(~exist('plxMat','var') || isempty(plxMat))
    [FileName,PathName] = uigetfile({'*.plx;*.PLX'},'Load plexon file');
    plxMat = fullfile(PathName,FileName);
end

if(~exist('NHP','var') || isempty(NHP))
    if(isstruct(plxMat))
        NHP = plxMat.File(1);
    else
        [~,cnm] = fileparts(plxMat);
        NHP = cnm(1); % assume default that first letter specifies NHP
    end
end

if(~exist('keep_sess','var') || isempty(keep_sess))
    keep_sess = 0; % [Hz] - keep times at this time resolution, but only use integer numbers for space reasons
end

if(~exist('skip_waves','var') || isempty(skip_waves))
    skip_waves = 0; % [Hz] - keep single spike waves
end

if(~exist('keep_nonstim','var') || isempty(keep_nonstim))
    keep_nonstim = 0; % [Hz] - keep single spike waves
end

% get event code information
if(~exist('EV','var') || isempty(EV))
    [FileName,PathName] = uigetfile({'*.pro;*.mat;*.m'},'Event code definition file');
    EV = fullfile(PathName,FileName);
end

if(~isstruct(EV))
     [~,~,cflext] = fileparts(EV);   
    
     if(strcmpi(cflext,'.pro'))
        [EV] = TEMPO_get_eventdefs(EV);
     elseif(strcmpi(cflext,'.mat'))
         EV = load(EV);
     elseif(strcmpi(cflext,'.m'))
         [evpth, evstem] = fileparts(EV);
         if(~exist(evstem,'file'))  % stupid matlab can not run functions unless they are in the known paths...
             addpath(evpth);
         end
         eval(sprintf('EV = %s;', evstem)); 
     else
         error('Event file format not recognized!');
     end
end

% ____________________________________________________________________________%
%% make date/animal dependent adjustments
if(NHP == 'E' || NHP == 'D')
    EV.TrialStart_ = EV.FixSpotOn_;
end

% ____________________________________________________________________________%
%% Get the Plexon data
if(~isstruct(plxMat))
    plxMat = PLX2Mat(plxMat, EV, keep_sess, skip_waves, PreTrialTime, PostTrialTime, 0, SampleRate);
end

%% determine stimulus onset
% plxMat.Task.StimOnsetToTrial = nan(plxMat.NTrials,1);
plxMat.Task.StimOnsetToTrial = plxMat.TrialMat_PDtimes(:,1);

% for(t=1:plxMat.NTrials)   
%     % PD event for stimulus onset
%     cpd = find(plxMat.TrialMat_PDtimes(t,:) >= plxMat.TrialStartEndDur(t,1), 1, 'first');
%     if(~isempty(cpd))
%         plxMat.Task.StimOnsetToTrial(t) = plxMat.TrialMat_PDtimes(t,cpd);
%     end
% end

%% remove trials without stimulus presentation
if(keep_nonstim == 0 && ~isempty(plxMat))
    rem_trial = [];
    for(t=1:plxMat.NTrials)
        p = find(plxMat.TrialMat_EVcodes(t,:) == EV.Target_, 1, 'first');
        if(isempty(p))
            rem_trial = [rem_trial; t]; %#ok<AGROW>
        end
    end
    
    plxMat.TrialNr = 1:plxMat.NTrials;
    
    if(~isempty(rem_trial))
        plxMat.NTrials = plxMat.NTrials - length(rem_trial);
        plxMat.TrialNr(rem_trial) = [];
        
        plxMat.TrialMat_INFOcodes(rem_trial,:) = [];
        plxMat.TrialMat_EVcodes(rem_trial,:) = [];
        plxMat.TrialMat_EVtimes(rem_trial,:) = [];
        plxMat.TrialMat_PDtimes(rem_trial,:) = [];
        
        plxMat.TrialStartEndDur(rem_trial,:) = [];
        
        plxMat.SessTime(rem_trial) = [];
        plxMat.Task.StimOnsetToTrial(rem_trial) = [];
        
        for(a = 1:length(plxMat.ActAD))
            cAD = sprintf('AD%02d', plxMat.ADChannel(a));
            if(isfield(plxMat,cAD))
                plxMat.(cAD)(rem_trial,:) = [];
            end
        end  

        for(u = 1:length(plxMat.DSP_names))
            if(isfield(plxMat.DSP,plxMat.DSP_names(u)))
                plxMat.DSP.(cell2mat(plxMat.DSP_names(u)))(rem_trial,:) = [];
            end
        end
    end
end

if(isempty(plxMat))
    warning('Not sufficient trials!');
    return;
end

% if(~exist('rig','var') || isempty(rig))
%     rig = plxMat.rig;
% end

% ____________________________________________________________________________%
%% get channel configuration

% this function refers to a script where changes specific to the set-up should be kept.
ADchan = PLX_define_ADchan(rig, NHP, plxMat.DateNum);
plxMat.NHP = NHP;

% ____________________________________________________________________________%
%% organize data according to event types and conditions.

% currently it is assumed that all data in the TrialMat_* is aligned to the stimulus
% onset as deside by the photo diode event.

% ____________________________________________________________________________%
%% get index vector for correct trials

% check if there are any abort codes and crop the info codes if needed
for(t=1:plxMat.NTrials)
    ep = find(plxMat.TrialMat_INFOcodes(t,:) == EV.Abort, 1, 'first');
    if(~isempty(ep) && ep<size(plxMat.TrialMat_INFOcodes,2))
        plxMat.TrialMat_INFOcodes(t,ep+1:end) = NaN;
    end
end

if(~any(isfinite(plxMat.TrialMat_INFOcodes(:,end))))
    plxMat.TrialMat_INFOcodes = wz_cropNaN(plxMat.TrialMat_INFOcodes);
end

% identify correct trials
plxMat.Task.Correct = nan(plxMat.NTrials,1);
for(t=1:plxMat.NTrials)
    evp = find(plxMat.TrialMat_INFOcodes(t,:) == EV.Correct_,1,'first');
    if(~isempty(evp))
        plxMat.Task.Correct(t) = plxMat.TrialMat_INFOcodes(t,evp+1);
    end
end

plxMat.NCorrect      = sum(plxMat.Task.Correct);
plxMat.Task.NTrials  = plxMat.NTrials;
plxMat.Task.NCorrect = plxMat.NCorrect;

% ____________________________________________________________________________%
%% Check for data format and task type

% Encodes for fixation spot onset/offset were later introduced (NHP D&E)

% are there any MG task specific encodes?
if(any(plxMat.TrialMat_INFOcodes(:) == EV.MG_Hold))
    MG_task = 1;
else
    MG_task = 0;
end

% are there any Search task specific encodes?
if(any(plxMat.TrialMat_INFOcodes(:) == EV.SetSize))
    Search_task = 1;
else
    Search_task = 0;
end

% are there any SAT task specific encodes?
if(any(plxMat.TrialMat_INFOcodes(:) == EV.SAT_Cond))
    SAT_task = 1;
else
    SAT_task = 0;
end

% Allocate trial information vectors
plxMat.Task.SRT              = nan(plxMat.NTrials,1);
plxMat.Task.SaccEnd          = nan(plxMat.NTrials,1);
plxMat.Task.SaccDur          = nan(plxMat.NTrials,1);
plxMat.Task.Reward           = nan(plxMat.NTrials,1);
plxMat.Task.Tone             = nan(plxMat.NTrials,1);
plxMat.Task.RewardTone       = nan(plxMat.NTrials,1);
plxMat.Task.ErrorTone        = nan(plxMat.NTrials,1);

if(any(plxMat.TrialMat_INFOcodes(:) == EV.FixSpotOn_))
    plxMat.Task.FixSpotOnTEMPO  = nan(plxMat.NTrials,1);
    plxMat.Task.FixSpotOffTEMPO = nan(plxMat.NTrials,1);
end

if(Search_task == 1)
    plxMat.Task.TargetLoc    = nan(plxMat.NTrials,1);
    plxMat.Task.TargetCol    = nan(plxMat.NTrials,1);
    plxMat.Task.IsCatch      = nan(plxMat.NTrials,1);

    plxMat.Task.TargetClut   = nan(plxMat.NTrials,1); % Min Target Color (CLUT value; higher values = BRIGHTER [more white])
    plxMat.Task.SetSize      = nan(plxMat.NTrials,1); % code meaning: 0 -> set size = 2; 1 -> set size = 4; 2 -> set size = 8
    plxMat.Task.SetSizeCond  = nan(plxMat.NTrials,1); % SET SIZE CONDITION (0 = 2; 1 = 4; 2 = 8; 3 = random (whatever this is supposed to mean)
    plxMat.Task.TaskType     = nan(plxMat.NTrials,1); % 0 = fixation; 1 = detection; 2 = search/MG??
    plxMat.Task.HoldTime     = nan(plxMat.NTrials,1); % use to determine search or memory guided (whatever this is supposed to mean)
    plxMat.Task.DistHomo     = nan(plxMat.NTrials,1); % Homogeneous (0)/Non-homogeneous (1)
    plxMat.Task.Eccentricity = nan(plxMat.NTrials,1);
end

plxMat.Task.FixTime      = nan(plxMat.NTrials,1); % FixTime_Jit_ is the ms value monkey must be holding fixation (with some uniform or exponentially distributed jitter) before the target appears
if(MG_task == 1)
    plxMat.Task.MG_Hold  = nan(plxMat.NTrials,1); % Code actual MG hold time
end

% Saccade Location as determined by TEMPO
if(~isempty(intersect(EV.SacLocVec, plxMat.TrialMat_INFOcodes(:))))
   plxMat.Task.SacLoc = nan(plxMat.NTrials,1);
end

% Stimuli Presented (needs to be checked what these codes are doing, maybe it could be simplified a lot)
if(~isempty(intersect(EV.StimLocVec, plxMat.TrialMat_INFOcodes(:))))
   plxMat.Task.StimTyp = nan(plxMat.NTrials,1);
end

% get some SAT specific information
if(SAT_task == 1)
    plxMat.Task.SAT_Block      = nan(plxMat.NTrials,1); % SAT blocking condtion (1 = blocked; 0 = random)
    plxMat.Task.SAT_Cond       = nan(plxMat.NTrials,1); % SAT condition (1 = slow, 2 = med, 3 = fast)
    plxMat.Task.SAT_CutOff     = nan(plxMat.NTrials,1); % SAT cutoff time for current block
    plxMat.Task.SAT_Percentile = nan(plxMat.NTrials,1); % Percentile (may by unnecessary)
    plxMat.Task.SAT_Reward     = nan(plxMat.NTrials,1); % Amount of reward for current block (time in ms solenoid on)
    plxMat.Task.SAT_Catch      = nan(plxMat.NTrials,1); % Catch Trial %
    plxMat.Task.SAT_PunishLoc  = nan(plxMat.NTrials,1); % Punish Time for Direction Error
    plxMat.Task.SAT_PunishTime = nan(plxMat.NTrials,1); % Punish Time for missed deadline (noResp)
    plxMat.Task.SAT_SearchDur  = nan(plxMat.NTrials,1); % Duration of SEARCH array (ms)
    plxMat.Task.SAT_FlipCond   = nan(plxMat.NTrials,1); % Flip Conditions automatically every X trials
    plxMat.Task.SAT_DispClear  = nan(plxMat.NTrials,1); % Was display cleared at deadline in FAST?
end

if(any(plxMat.TrialMat_INFOcodes(:) == EV.Stimulation_))
    Task.StimEV    = nan(plxMat.NTrials,1);
    Task.StimTrial = nan(plxMat.NTrials,1);
    Task.StimProb  = nan(plxMat.NTrials,1);
    Task.StimDur   = nan(plxMat.NTrials,1);
    Task.StimTime  = nan(plxMat.NTrials,1);   
end  % for(t=1:plxMat.NTrials)

plxMat.Task.Gains_EyeX    = nan(plxMat.NTrials,1);
plxMat.Task.Gains_EyeY    = nan(plxMat.NTrials,1);
plxMat.Task.Gains_XYF     = nan(plxMat.NTrials,1);
plxMat.Task.Gains_YXF     = nan(plxMat.NTrials,1);

% loop over trials and identify relevant information
for(t=1:plxMat.NTrials)
    cInf = plxMat.TrialMat_INFOcodes(t,:);  % get the info codes corresponding to the current trial

    % in this data format the event encodes are used in very inconsistent
    % ways. One variant uses a event marker to identify the information type and
    % the subsequent number specifies the value for this information. The
    % second variant uses a arbitrary number with the information number
    % added to it. This inconsistent coding makes it difficult to come up
    % with a clean and straightforward interpretation method.
    % also, some timing informations are kept as encoded numbers as
    % determined by tempo, not as actual event codes at the corresponding
    % times. This is a dangerous method because there is no way to validate
    % these numbers or align them with other information such as screen
    % refresh. My humble recommendation, handle this data, if you can not avoid it,
    % with a lot o caution.

    plxMat.Task.TargetLoc(t)    = PLXin_get_event_value(EV.TargetLoc, t);
    plxMat.Task.TargetCol(t)    = PLXin_get_event_value(EV.TargetLoc, t);
    plxMat.Task.IsCatch(t)      = plxMat.Task.TargetLoc(t) == 255;
     
    plxMat.Task.TargetClut(t)   = PLXin_get_event_value(EV.TargetClut,   t); % Min Target Color (CLUT value; higher values = BRIGHTER [more white])
    plxMat.Task.SetSize(t)      = PLXin_get_event_value(EV.SetSize,      t); % code meaning: 0 -> set size = 2; 1 -> set size = 4; 2 -> set size = 8
    plxMat.Task.SetSizeCond(t)  = PLXin_get_event_value(EV.SetSizeCond,  t); % SET SIZE CONDITION (0 = 2; 1 = 4; 2 = 8; 3 = random (whatever this is supposed to mean)
    plxMat.Task.TaskType(t)     = PLXin_get_event_value(EV.TaskType,     t); % 0 = fixation; 1 = detection; 2 = search/MG??
    plxMat.Task.HoldTime(t)     = PLXin_get_event_value(EV.HoldTime,     t); % use to determine search or memory guided (whatever this is supposed to mean)
    plxMat.Task.DistHomo(t)     = PLXin_get_event_value(EV.DistHomo,     t); % Homogeneous (0)/Non-homogeneous (1)
    plxMat.Task.Eccentricity(t) = PLXin_get_event_value(EV.Eccentricity, t);
   
    plxMat.Task.FixTime(t)     = PLXin_get_event_value(EV.FixTime,       t); % FixTime_Jit_ is the ms value monkey must be holding fixation (with some uniform or exponentially distributed jitter) before the target appears
    if(MG_task == 1)
        plxMat.Task.MG_Hold(t) = PLXin_get_event_value(EV.MG_Hold,       t); % Code actual MG hold time
    end

   % Saccade Location as determined by TEMPO
   [~, sp] = intersect(EV.SacLocVec, cInf);
   if(~isempty(sp))
       if(length(sp) > 1) 
           warning('Something strange happened, more than one saccade target encoded!');
           sp = sp(1);
       end
       plxMat.Task.SacLoc(t) = sp-1;
   end

    % Stimuli Presented (needs to be checked what these codes are doing, maybe it could be simplified a lot)
    [~, sl, el] = intersect(EV.StimLocVec, cInf);
    if(~isempty(el))
       plxMat.Task.StimTyp(t,sl) = cInf(el+1);
    end

    % get some SAT specific information
    if(SAT_task == 1)
        plxMat.Task.SAT_Block(t)      = PLXin_get_event_value(EV.SAT_Block,      t); % SAT blocking condtion (1 = blocked; 0 = random)
        plxMat.Task.SAT_Cond(t)       = PLXin_get_event_value(EV.SAT_Cond,       t); % SAT condition (1 = slow, 2 = med, 3 = fast)
        plxMat.Task.SAT_CutOff(t)     = PLXin_get_event_value(EV.SAT_CutOff,     t); % SAT cutoff time for current block
        plxMat.Task.SAT_Percentile(t) = PLXin_get_event_value(EV.SAT_Percentile, t); % Percentile (may by unnecessary)
        plxMat.Task.SAT_Reward(t)     = PLXin_get_event_value(EV.SAT_Reward,     t); % Amount of reward for current block (time in ms solenoid on)
        plxMat.Task.SAT_Catch(t)      = PLXin_get_event_value(EV.SAT_Catch,      t); % Catch Trial %
        plxMat.Task.SAT_PunishLoc(t)  = PLXin_get_event_value(EV.SAT_PunishLoc,  t); % Punish Time for Direction Error
        plxMat.Task.SAT_PunishTime(t) = PLXin_get_event_value(EV.SAT_PunishTime, t); % Punish Time for missed deadline (noResp)
        plxMat.Task.SAT_SearchDur(t)  = PLXin_get_event_value(EV.SAT_SearchDur,  t); % Duration of SEARCH array (ms)
        plxMat.Task.SAT_FlipCond(t)   = PLXin_get_event_value(EV.SAT_FlipCond,   t); % Flip Conditions automatically every X trials
        plxMat.Task.SAT_DispClear(t)  = PLXin_get_event_time(EV.SAT_DispClear,   t); % Was display cleared at deadline in FAST?
    end
    
    % eye channel information
    plxMat.Task.Gains_EyeX(t) = CorrGain(PLXin_get_event_time(EV.GainsX,    t));
    plxMat.Task.Gains_EyeY(t) = CorrGain(PLXin_get_event_time(EV.GainsY,    t));
    plxMat.Task.Gains_XYF(t)  = CorrGain(PLXin_get_event_time(EV.Gains_XYF, t));
    plxMat.Task.Gains_YXF(t)  = CorrGain(PLXin_get_event_time(EV.Gains_YXF, t));
   
    % check for the ocurrence of stimulation       
    Task.StimEV    = nan(plxMat.NTrials,1);
    Task.StimProb  = nan(plxMat.NTrials,1);
    Task.StimDur   = nan(plxMat.NTrials,1);
    Task.StimTime  = nan(plxMat.NTrials,1);   
    
    if(isfield(plxMat.Task, 'StimEV'))
        stimpos = find(plxMat.TrialMat_EVcodes(t,:) == 666);
        if(~isempty(stimpos))  
            plxMat.Task.StimTrial(t) = 1;
            if(length(stimpos) > 1)
                warning(['More than one stimulation occurrence in trial ',int2str(t),'! Only log the first one and ignore the rest!']);
                stimpos(2:end) = [];
            end
            plxMat.Task.StimEV(t) = plxMat.TrialMat_EVtimes(t,stimpos); 
        else
            % no stimulation in this trial
            plxMat.Task.StimTrial(t) = 0;
        end
    end
end  % for(t=1:plxMat.NTrials)

% correct for cases where the abort code just occurs after the FixTime code
plxMat.Task.FixTime(plxMat.Task.FixTime == EV.Abort) = NaN;

% convert into logical
if(SAT_task == 1)
    plxMat.Task.SAT_DispClear = isfinite(plxMat.Task.SAT_DispClear);
end

% ____________________________________________________________________________%
%% define error outcome
% This needs to be checked and coded more dynamically
plxMat.Task.error_names = {'False', 'Early', 'Late', 'FixBreak', 'HoldError', ...
                           'CatchError', 'Correct_too_fast', 'Correct_too_slow'};

plxMat.Task.error = nan(plxMat.NTrials,1);
plxMat.Task.error(plxMat.Task.Correct==1) = 0;

false_resp = logical(sum(bsxfun(@eq, plxMat.TrialMat_EVcodes, EV.Error_sacc), 2));
if(any(false_resp))
    plxMat.Task.error(false_resp) = 1;
end

early_resp = logical(sum(bsxfun(@eq, plxMat.TrialMat_EVcodes, EV.HoldError), 2));
if(any(early_resp))
    plxMat.Task.error(early_resp) = 2;
end

fix_break = logical(sum(bsxfun(@eq, plxMat.TrialMat_EVcodes, EV.FixError_), 2));
if(any(fix_break))
    plxMat.Task.error(fix_break) = 4;
end

hold_err = logical(sum(bsxfun(@eq, plxMat.TrialMat_EVcodes, EV.TargetHold_Error), 2));
if(any(hold_err))
    plxMat.Task.error(hold_err) = 5;
end

catch_go = logical(sum(bsxfun(@eq, plxMat.TrialMat_EVcodes, EV.CatchError), 2));
if(any(catch_go))
    plxMat.Task.error(catch_go) = 6;
end

catch_hold = logical(sum(bsxfun(@eq, plxMat.TrialMat_EVcodes, EV.Correct_too_fast), 2));
if(any(catch_hold))
    plxMat.Task.error(catch_hold) = 7;
end

catch_hold = logical(sum(bsxfun(@eq, plxMat.TrialMat_EVcodes, EV.Correct_too_slow), 2));
if(any(catch_hold))
    plxMat.Task.error(catch_hold) = 8;
end

% if no error codes were found try to determine trial error type
if(any(plxMat.Task.error ~= 0 & isfinite(plxMat.Task.error) == 1) == 0)
    warning('Apparently no error codes in data file. Trying to reconstruct');
    notcorrect = find(plxMat.Task.Correct == 0);
    
    for(et = 1:length(notcorrect))
        cet = notcorrect(et);
        % saccade made, item hold sufficiently long, but wrong item selected
        if(plxMat.Task.SRT(cet) >  0 && any(plxMat.TrialMat_EVcodes(cet,:) == EV.Decide_))
            plxMat.Task.error(cet) = 1;
        % saccade made but item not fixated sufficiently long (unclear if correct or incorrect choice)
        elseif(plxMat.Task.SRT(cet) >  0 ) 
            plxMat.Task.error(cet) = 5;
        % fixation break prior to stimulus onset    
        elseif(isnan(plxMat.Task.StimOnsetToTrial(t)))
            plxMat.Task.error(cet) = 4;
        % fixation break after stimulus onset but before go signal
        elseif(any(plxMat.TrialMat_EVcodes(cet,:) == EV.Fixate_) && isfinite(plxMat.Task.StimOnsetToTrial(cet)) && ...
               isnan(plxMat.Task.GoCue(cet)) && isnan(plxMat.Task.SRT(cet))) 
            plxMat.Task.error(cet) = 2;
        % no sacade made in time after go signal
        elseif(any(plxMat.TrialMat_EVcodes(cet,:) == EV.Fixate_) && isfinite(plxMat.Task.GoCue(cet)) && ...
               isnan(plxMat.Task.SRT(cet))) 
            plxMat.Task.error(cet) = 3;
        % trial abort due to fixation break    
        elseif(any(plxMat.TrialMat_EVcodes(cet,:) == EV.Abort))
            plxMat.Task.error(cet) = 4;
        elseif(any(plxMat.TrialMat_EVcodes(cet,:) == EV.Abort_old))
            plxMat.Task.error(cet) = 4;
        else
            warning('Not able to identify error type');
        end
    end
end

% ____________________________________________________________________________%
%% Clean up analog channels and spike data

% align data now to stimulus onset.
timebin_mat    = bsxfun(@minus, plxMat.timevec, round(plxMat.Task.StimOnsetToTrial));
plxMat.timevec = PreOnsetTime : max(timebin_mat(:));
AD_mat = nan(plxMat.NTrials, length(plxMat.timevec)); % use to pre-allocate data

% ____________________________________________________________________________%
%% clean spike waveforms if they should not be kept
% if(skip_waves == 1) % remove single wave forms if they are still in the data and not wanted.
%     if(isfield(plxMat.DSP.([cell2mat(plxMat.DSP_names(1)),'_wave']),'waves'))
%         for(s=1:length(plxMat.DSP_names))
%             if(isfield(plxMat.DSP.([cell2mat(plxMat.DSP_names(s)),'_wave']),'waves'))
%                 plxMat.DSP.([cell2mat(plxMat.DSP_names(s)),'_wave']) = ...
%                     rmfield(plxMat.DSP.([cell2mat(plxMat.DSP_names(s)),'_wave']),'waves');
%             end
%         end
%     end
% end

% ____________________________________________________________________________%
% remove irrelevant spike times
for(s=1:length(plxMat.DSP_names))
    cspktimes = bsxfun(@minus, plxMat.DSP.(cell2mat(plxMat.DSP_names(s))), plxMat.Task.StimOnsetToTrial);
    cspktimes(cspktimes < PreOnsetTime | cspktimes > plxMat.timevec(end)) = NaN;
    plxMat.DSP.(cell2mat(plxMat.DSP_names(s))) = wz_cropNaN(cspktimes);
end

% ____________________________________________________________________________%
%% identify LFP channels
switch lower(plxMat.ElectrodeType)
    case 'neuronexus';
        SpikeChan = 1:32;
    case 'uprobe'
        SpikeChan = 1:24;
    case {'single', 'array', 'na'}
        SpikeChan = plxMat.SpikeChan; % better define it extra in the ADchan configuration to allow LFP data without spikes.
end

plxMat.LFPchan = SpikeChan; % this needs some adjustment to work reliable.
plxMat.num_lfp = length(plxMat.LFPchan);
plxMat.LFP_names = {};

for(s=1:length(SpikeChan))
    spk_chan = SpikeChan(s);
    cAD  = sprintf('AD%02d',  spk_chan);
    cLFP = sprintf('LFP%02d', spk_chan);

    plxMat.LFP_names(s,1) = {cLFP};

    gotLFP = 0;
    
    if(isfield(plxMat,cAD))
        gotLFP = 1;
    elseif(isfield(plxMat,'LFP'))
        if(isfield(plxMat.LFP,cLFP))
            gotLFP = 1;
        end
    end
        
    if(gotLFP == 0)
        warning('Did not find a corresponding analog channel for the current spike channel!');
    else
        % pre-allocate
        plxMat.LFP.(cLFP) = AD_mat;
        
        % just another painful trial loop...
        for(ct = 1:plxMat.NTrials)
            cpos = find(timebin_mat(ct,:) >= PreOnsetTime);
            [~, spos] = intersect(plxMat.timevec, timebin_mat(ct,cpos));
            
            if(~isempty(cpos)) % discard data from trials without stimulus onset
                plxMat.LFP.(cLFP)(ct,spos) = plxMat.(cAD)(ct,cpos);
            end
        end
        plxMat   =  rmfield(plxMat, cAD);
%         plxMat.LFP.([cLFP,'_z']) = plxMat.([cAD,'_z']);
%         plxMat   =  rmfield(plxMat, [cAD,'_z']);

        chanpos = find(plxMat.ADChannel == spk_chan - 1);
        plxMat.LFP.([cLFP,'_info']).orig      = plxMat.ADName(    chanpos);
        plxMat.LFP.([cLFP,'_info']).gain      = plxMat.ADGain(    chanpos);
        plxMat.LFP.([cLFP,'_info']).pregain   = plxMat.ADPreGain( chanpos);
        plxMat.LFP.([cLFP,'_info']).channel   = plxMat.ADChannel( chanpos);
        plxMat.LFP.([cLFP,'_info']).rate      = plxMat.ADSmplRate(chanpos);
        plxMat.LFP.([cLFP,'_info']).timestamp = plxMat.ADStartRec(chanpos); % this should be corrected by now

        if(isfield(plxMat,[cAD,'_session']))
            if(keep_sess == 1)
                plxMat.LFP.([cLFP,'_session']) = plxMat.([cAD,'_session']);
            end
            plxMat = rmfield(plxMat, [cAD,'_session']);
        end
    end
end

% ____________________________________________________________________________%
%% get Eye data
plxMat.Eye.EyeX  = [];
plxMat.Eye.EyeY  = [];
plxMat.Eye.Pupil = [];

% x-coordinate of eye position
if(~isempty(ADchan.EyeX) && ~isempty(ADchan.EyeY))
    eyeX_chan = ADchan.EyeX;
    eX = sprintf('AD%02d', eyeX_chan);
    
    eyeY_chan = ADchan.EyeY;
    eY = sprintf('AD%02d', eyeY_chan);
    
    if(~isfield(plxMat, eX) && ~isfield(plxMat.Eye, 'EyeX'))
        warning('Did not find a corresponding analog channel for the X eye position!');
    else
        % pre-allocate
        plxMat.Eye.EyeX = AD_mat;
        plxMat.Eye.EyeY = AD_mat;
        
        % just another painful trial loop...
        for(ct = 1:plxMat.NTrials)
            cpos = find(timebin_mat(ct,:) >= PreOnsetTime);
            [~, spos] = intersect(plxMat.timevec, timebin_mat(ct,cpos));
            
            if(~isempty(cpos)) % discard data from trials without stimulus onset
                plxMat.Eye.EyeX(ct,spos) = plxMat.(eX)(ct,cpos);
                plxMat.Eye.EyeY(ct,spos) = plxMat.(eY)(ct,cpos);
                
            end
        end
        plxMat = rmfield(plxMat, eX);
        plxMat = rmfield(plxMat, eY);
        
        chanpos = find(plxMat.ADChannel == eyeX_chan - 1);
        plxMat.Eye.EyeX_info.orig      = plxMat.ADName(    chanpos);
        plxMat.Eye.EyeX_info.gain      = plxMat.ADGain(    chanpos);
        plxMat.Eye.EyeX_info.pregain   = plxMat.ADPreGain( chanpos);
        plxMat.Eye.EyeX_info.channel   = plxMat.ADChannel( chanpos);
        plxMat.Eye.EyeX_info.rate      = plxMat.ADSmplRate(chanpos);
        plxMat.Eye.EyeX_info.timestamp = plxMat.ADStartRec(chanpos); % this should be corrected by now
        
        chanpos = find(plxMat.ADChannel == eyeY_chan - 1);
        plxMat.Eye.EyeY_info.orig      = plxMat.ADName(    chanpos);
        plxMat.Eye.EyeY_info.gain      = plxMat.ADGain(    chanpos);
        plxMat.Eye.EyeY_info.pregain   = plxMat.ADPreGain( chanpos);
        plxMat.Eye.EyeY_info.channel   = plxMat.ADChannel( chanpos);
        plxMat.Eye.EyeY_info.rate      = plxMat.ADSmplRate(chanpos);
        plxMat.Eye.EyeY_info.timestamp = plxMat.ADStartRec(chanpos); % this should be corrected by now
        
    end
end

% calculate distance of subsequent gaze positions
if(~isempty(plxMat.Eye.EyeX) && ~isempty(plxMat.Eye.EyeY))
    plxMat.Eye.Dist = nan(size(plxMat.Eye.EyeY));
    plxMat.Eye.Dist(:,2:end) = sqrt(diff(plxMat.Eye.EyeX,1,2).^2 + diff(plxMat.Eye.EyeY,1,2).^2);
end

% Pupil size
if(~isempty(ADchan.Pupil))
    eye_chan = ADchan.Pupil;
    eP = sprintf('AD%02d', eye_chan);
    
    if(~isfield(plxMat, eP) && ~isfield(plxMat.Eye, 'Pupil'))
        warning('Did not find a corresponding analog channel for the pupil size!');
    else
        % pre-allocate
        plxMat.Eye.Pupil = AD_mat;
        
        % just another painful trial loop...
        for(ct = 1:plxMat.NTrials)
            cpos = find(timebin_mat(ct,:) >= PreOnsetTime);
            [~, spos] = intersect(plxMat.timevec, timebin_mat(ct,cpos));
            
            if(~isempty(cpos)) % discard data from trials without stimulus onset
                plxMat.Eye.Pupil(ct,spos) = plxMat.(eP)(ct,cpos);
            end
        end
        plxMat = rmfield(plxMat, eP);
        %         plxMat.EEG.Pupil_z = plxMat.([eP,'_z']);
        %         plxMat = rmfield(plxMat, [eP,'_z']);
        
        chanpos = find(plxMat.ADChannel == eye_chan - 1);
        plxMat.Eye.Pupil_info.orig      = plxMat.ADName(    chanpos);
        plxMat.Eye.Pupil_info.gain      = plxMat.ADGain(    chanpos);
        plxMat.Eye.Pupil_info.pregain   = plxMat.ADPreGain( chanpos);
        plxMat.Eye.Pupil_info.channel   = plxMat.ADChannel( chanpos);
        plxMat.Eye.Pupil_info.rate      = plxMat.ADSmplRate(chanpos);
        plxMat.Eye.Pupil_info.timestamp = plxMat.ADStartRec(chanpos); % this should be corrected by now
    end
end

% ____________________________________________________________________________%
%% Determine data size
if(exist('ByetSize','file'))
    fprintf('Resulting structure has a size of %s', ByetSize(plxMat));
end

% ____________________________________________________________________________%
%% helper functions

function ctm = PLXin_get_event_time(EVcode, ctrial)
% look for the occurrence of the specified event code and return the
% corresponding event time
    % use of evalin usually not recommended. Here, as kind of inline code
    EVc = evalin('caller',['plxMat.TrialMat_EVcodes(',int2str(ctrial),',:)']);
    EVt = evalin('caller',['plxMat.TrialMat_EVtimes(',int2str(ctrial),',:)']);
    ctm = NaN;

    p = EVc == EVcode;
    if(sum(p)==1)
        ctm = EVt(p);
    elseif(sum(p)>1)
        warning('Something weird occurred: more than one event found (only first will be used). Check carefully!');
        ctm = ctm(1);
    end

function cval = PLXin_get_event_value(EVcode, ctrial, cpos)
% look for the occurrence of the specified event code and return the
% value of one of the subsequent event code position
    if(~exist('cpos','var') || isempty(cpos))
        cpos = 1; 
    end

    % use of evalin usually not recommended. Here, as kind of inline code
    EVc = evalin('caller',['plxMat.TrialMat_EVcodes(',int2str(ctrial),',:)']);
    EVt = evalin('caller',['plxMat.TrialMat_EVtimes(',int2str(ctrial),',:)']);
    cval = NaN;

    p = find(EVc == EVcode);
    if(length(p) > 1)
        warning('Something weird occurred: more than one event found. Check carefully!');
        p = p(1);
    end

    if(~isempty(p) && p+cpos <= length(EVc))
        cval = EVt(p+cpos);
    end

    
function  corval = CorrGain(val)
% apparently the definition of gain values did change in the data format.
% This function reads out the encoded value and applies the hopefully
% appropriate correction.

    % val = PLXin_get_event_time(EVcode, ctrial);

    if(val > 30000)
        corval = -1 * (val - 30000);
    else
        corval =      (val - 20000);        
    end
