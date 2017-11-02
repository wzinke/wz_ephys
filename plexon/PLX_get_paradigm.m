function plxMat = PLX_get_paradigm(plxMat, EV, rig, NHP, keep_sess, skip_waves, keep_nonstim, ITIstim)
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
%       ITIstim       - check for stimulation codes during the inter trial interval
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
PreTrialTime  =  -50; % start sampling trial based data at this time relative to trial onset.
PostTrialTime =  500; % continue to collect analog data for this time period - keep it longer to keep reward related response modulations.
PreOnsetTime  = -500; % discard analog data before this time, pad with nans if this time is not matched
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

if(~exist('ITIstim','var') || isempty(ITIstim))
    ITIstim = 0; % flag to indicate whether stimulation occurred between trial (ad hoc hack for ultrasound experiments)
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
%% Get the Plexon data
if(~isstruct(plxMat))
    plxMat = PLX2Mat(plxMat, EV, keep_sess, skip_waves, PreTrialTime, PostTrialTime, ITIstim, SampleRate);
end

%% determine stimulus onset
plxMat.Task.StimOnsetToTrial = nan(plxMat.NTrials,1);

for(t=1:plxMat.NTrials)
    tempoStimOn = PLXin_get_event_time(EV.Target_, t);

    % PD event for stimulus onset
    cpd = find(plxMat.TrialMat_PDtimes(t,:) >= tempoStimOn, 1, 'first');
    if(~isempty(cpd))
        plxMat.Task.StimOnsetToTrial(t) = plxMat.TrialMat_PDtimes(t,cpd);
    end
end

%% remove trials without stimulus presentation
if(keep_nonstim == 0 && ~isempty(plxMat))
    rem_trial = [];
    for(t=1:plxMat.NTrials)
        p = find(plxMat.TrialMat_EVcodes(t,:) == EV.Target_, 1, 'first');
        if(isempty(p))
            rem_trial = [rem_trial; t];
        end
    end

    plxMat.TrialNr = 1:plxMat.NTrials;

    if(~isempty(rem_trial))
        plxMat.NTrials = plxMat.NTrials - length(rem_trial);
        plxMat.TrialNr(rem_trial) = [];

        plxMat.TrialMat_INFOcodes(rem_trial,:)  = [];
        plxMat.TrialMat_EVcodes(rem_trial,:)    = [];
        plxMat.TrialMat_EVtimes(rem_trial,:)    = [];
        plxMat.TrialMat_PDtimes(rem_trial,:)    = [];

        plxMat.TrialStartEndDur(rem_trial,:)    = [];

        plxMat.Task.StimOnsetToTrial(rem_trial) = [];

        plxMat.SessTime(rem_trial) = [];

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
plxMat.Task.Correct  = logical(sum(bsxfun(@eq, plxMat.TrialMat_EVcodes, EV.Correct_),2));
plxMat.NCorrect      = sum(plxMat.Task.Correct);
plxMat.Task.NTrials  = plxMat.NTrials;
plxMat.Task.NCorrect = plxMat.NCorrect;

% ____________________________________________________________________________%
%% Stimulus onset time (relative to trial start, i.e. first time bin of analog data

% Allocate trial information vectors
plxMat.Task.SRT              = nan(plxMat.NTrials,1);  % saccadic response time relative to go signal
plxMat.Task.Saccade          = nan(plxMat.NTrials,1);  % time of saccade relative to trial start
plxMat.Task.SaccEnd          = nan(plxMat.NTrials,1);
plxMat.Task.SaccDur          = nan(plxMat.NTrials,1);
plxMat.Task.Reward           = nan(plxMat.NTrials,1);
plxMat.Task.Tone             = nan(plxMat.NTrials,1);
plxMat.Task.RewardTone       = nan(plxMat.NTrials,1);
plxMat.Task.ErrorTone        = nan(plxMat.NTrials,1);
plxMat.Task.FixSpotOn        = nan(plxMat.NTrials,1);
plxMat.Task.FixSpotOff       = nan(plxMat.NTrials,1);
plxMat.Task.FixSpotOnTEMPO   = nan(plxMat.NTrials,1);
plxMat.Task.FixSpotOffTEMPO  = nan(plxMat.NTrials,1);

% ____________________________________________________________________________%
%% loop over trials to get relevant timings.

% Get photo diode event as stimulus onset, keep this time relative to trial
% onset for later time adjustments.

for(t=1:plxMat.NTrials)
    tempoStimOn = PLXin_get_event_time(EV.Target_, t);

    plxMat.Task.refresh_offset(t)  = plxMat.Task.StimOnsetToTrial(t) - tempoStimOn;
    plxMat.Task.Saccade(t)         = PLXin_get_event_time(EV.Saccade_,    t);
    plxMat.Task.SaccEnd(t)         = PLXin_get_event_time(EV.Decide_,     t);
    plxMat.Task.Reward(t)          = PLXin_get_event_time(EV.Reward_,     t);
    plxMat.Task.Tone(t)            = PLXin_get_event_time(EV.Tone_,       t);
    plxMat.Task.RewardTone(t)      = PLXin_get_event_time(EV.Reward_tone, t);
    plxMat.Task.ErrorTone(t)       = PLXin_get_event_time(EV.Error_tone,  t);
    plxMat.Task.FixSpotOnTEMPO(t)  = PLXin_get_event_time(EV.FixSpotOn_,  t);
    plxMat.Task.FixSpotOffTEMPO(t) = PLXin_get_event_time(EV.FixSpotOff_, t);

    % PD event for fixation spot offset (not working for MG task,
    % correction below)
    cpd = find(plxMat.TrialMat_PDtimes(t,:) >= plxMat.Task.FixSpotOffTEMPO(t), 1, 'first');
    if(~isempty(cpd))
        plxMat.Task.FixSpotOff(t) = plxMat.TrialMat_PDtimes(t,cpd);
    end

    % PD code for fixation spot onset
    cpd = find(plxMat.TrialMat_PDtimes(t,:) >= plxMat.Task.FixSpotOnTEMPO(t), 1, 'first');
    if(~isempty(cpd))
        plxMat.Task.FixSpotOn(t) = plxMat.TrialMat_PDtimes(t,cpd);
    end
end

% define stimulus onset as zero
%plxMat.Task.Task.SRT       = plxMat.Task.SRT              - plxMat.Task.StimOnsetToTrial;
plxMat.Task.Saccade         = plxMat.Task.Saccade          - plxMat.Task.StimOnsetToTrial;  % time of saccade relative to trial start
plxMat.Task.SaccEnd         = plxMat.Task.SaccEnd          - plxMat.Task.StimOnsetToTrial;
plxMat.Task.Reward          = plxMat.Task.Reward           - plxMat.Task.StimOnsetToTrial;
plxMat.Task.Tone            = plxMat.Task.Tone             - plxMat.Task.StimOnsetToTrial;
plxMat.Task.RewardTone      = plxMat.Task.RewardTone       - plxMat.Task.StimOnsetToTrial;
plxMat.Task.ErrorTone       = plxMat.Task.ErrorTone        - plxMat.Task.StimOnsetToTrial;
plxMat.Task.FixSpotOn       = plxMat.Task.FixSpotOn        - plxMat.Task.StimOnsetToTrial;
plxMat.Task.FixSpotOff      = plxMat.Task.FixSpotOff       - plxMat.Task.StimOnsetToTrial;
plxMat.Task.StimOnset       = plxMat.Task.StimOnsetToTrial - plxMat.Task.StimOnsetToTrial; % should be all zero aferwards
plxMat.Task.FixSpotOnTEMPO  = plxMat.Task.FixSpotOnTEMPO   - plxMat.Task.StimOnsetToTrial;
plxMat.Task.FixSpotOffTEMPO = plxMat.Task.FixSpotOffTEMPO  - plxMat.Task.StimOnsetToTrial;

plxMat.Task.SaccDur         = plxMat.Task.SaccEnd          - plxMat.Task.Saccade;


% keep event times in session time as well
if(keep_sess == 1)
%   plxMat.Task.Sess.SRT        = plxMat.TrialMat_EVtimes(plxMat.TrialMat_EVcodes == EV.Saccade_);
    plxMat.Task.Sess.Saccade    = plxMat.TrialMat_EVtimes(plxMat.TrialMat_EVcodes == EV.Saccade_);
    plxMat.Task.Sess.SaccEnd    = plxMat.TrialMat_EVtimes(plxMat.TrialMat_EVcodes == EV.Decide_);
    plxMat.Task.Sess.Reward     = plxMat.TrialMat_EVtimes(plxMat.TrialMat_EVcodes == EV.Reward_);
    plxMat.Task.Sess.Tone       = plxMat.TrialMat_EVtimes(plxMat.TrialMat_EVcodes == EV.Tone_);
    plxMat.Task.Sess.RewardTone = plxMat.TrialMat_EVtimes(plxMat.TrialMat_EVcodes == EV.Reward_tone);
    plxMat.Task.Sess.ErrorTone  = plxMat.TrialMat_EVtimes(plxMat.TrialMat_EVcodes == EV.Error_tone);
    plxMat.Task.Sess.FixSpotOn  = plxMat.TrialMat_EVtimes(plxMat.TrialMat_EVcodes == EV.FixSpotOn_);
%   plxMat.Task.Sess.FixSpotOff = plxMat.TrialMat_EVtimes(plxMat.TrialMat_EVcodes == EV.FixSpotOff_);
end

% ____________________________________________________________________________%
% loop over trials to get stimulus information
%% stimulus information

% this needs a lot of re-hacking. Ideally, some paradigm information will be
% included that defines code and parameters, maybe a kind of plug-in inline function
% that takes the event codes and stimulus information plot and returns condition
% labels and stimulus information.

% ____________________________________________________________________________%
%% Get trial information
% Here, the plxMat.TrialMat_INFOcodes matrix is interpreted. This matrix corresponds
% to the trial info block created by Tempo. Thus, the columns represent the
% different entries for the task parameters specified for Tempo. Check the
% INFOS.pro file of the tempo configuration for the interpretation of the columns.

% predefine expected Task variables
plxMat.Task.TargetLoc    = nan(plxMat.NTrials,1);  % Target location (angle)
plxMat.Task.TargetPos    = nan(plxMat.NTrials,1);  % Target location (index number)
plxMat.Task.Eccentricity = nan(plxMat.NTrials,1);  % Eccentricity of search array items
plxMat.Task.SOA          = nan(plxMat.NTrials,1);  % ?
plxMat.Task.TargetHold   = nan(plxMat.NTrials,1);  % time to keep gaze at target to be accepted as correct response
plxMat.Task.TrialType    = nan(plxMat.NTrials,1);  % Trial Type (random vs repeated displays)

plxMat.Task.GoCue         = nan(plxMat.NTrials,1);  % event that triggers the animal to respond

plxMat.Task.MaxSaccDur   = nan(plxMat.NTrials,1);
plxMat.Task.MaxSaccTime  = nan(plxMat.NTrials,1);
plxMat.Task.TimeOut      = nan(plxMat.NTrials,1);
plxMat.Task.RewardDur    = nan(plxMat.NTrials,1);
plxMat.Task.RewardOffset = nan(plxMat.NTrials,1);
plxMat.Task.TargetHold   = nan(plxMat.NTrials,1);
plxMat.Task.Gains_EyeX   = nan(plxMat.NTrials,1);
plxMat.Task.Gains_EyeY   = nan(plxMat.NTrials,1);
plxMat.Task.Gains_XYF    = nan(plxMat.NTrials,1);
plxMat.Task.Gains_YXF    = nan(plxMat.NTrials,1);

if(any(plxMat.TrialMat_INFOcodes(:,1) == EV.InfosZero))  % MG trials exist
    plxMat.Task.MGcolor       = nan(plxMat.NTrials,1);  % color of item in the MG task
    plxMat.Task.DelayDur      = nan(plxMat.NTrials,1);  % Delay duration for memory task
    plxMat.Task.TargetCol_R   = nan(plxMat.NTrials,1);  % RGB color of target item - R value
    plxMat.Task.TargetCol_G   = nan(plxMat.NTrials,1);  % RGB color of target item - G value
    plxMat.Task.TargetCol_B   = nan(plxMat.NTrials,1);  % RGB color of target item - B value
end

if(any(plxMat.TrialMat_INFOcodes(:) == 3999))  % Search task (including capture) trials exist
    plxMat.Task.ArrayStruct   = nan(plxMat.NTrials,1);  % Structured vs unstrcutured display (typical vs. contextual cue)
    plxMat.Task.DistHomo      = nan(plxMat.NTrials,1);  % homo/hetero/random search cond
    plxMat.Task.SingMode      = nan(plxMat.NTrials,1);  % Could a singleton be expected in a trial?
    plxMat.Task.SetSize       = nan(plxMat.NTrials,1);  % Number of items in search array
    plxMat.Task.TargetType    = nan(plxMat.NTrials,1);  % Target Type - T (1) vs L (2)
    plxMat.Task.Singleton     = nan(plxMat.NTrials,1);  % singleton presence: 1111 - singleton absent; 2222 - singleton present
    plxMat.Task.TargetHemi    = repmat('-',plxMat.NTrials,1);  % Target hemifield: 8200 - right; 8100 - left;8888 - midline
    plxMat.Task.DistHemi      = repmat('-',plxMat.NTrials,1);  % Distractor hemifield: 8200 - right; 8100 - left;8888 - midline
    plxMat.Task.DistLoc       = nan(plxMat.NTrials,1);  % Distractor location (angle)
    plxMat.Task.DistPos       = nan(plxMat.NTrials,1);  % Distractor location (index number)
    plxMat.Task.IsCatch       = nan(plxMat.NTrials,1);  % current trial is catch trial
    plxMat.Task.SingleDistCol = nan(plxMat.NTrials,1);  % Singleton Distractor Color
    plxMat.Task.DistOri       = nan(plxMat.NTrials,1);  % Distractor Orientation: 0-3 (0 = upright, 1= right tilt, 3 = inverted, 4 = left tilt)
    plxMat.Task.TargetOri     = nan(plxMat.NTrials,1);  % Target Orientation:     0-3 (0 = upright, 1= right tilt, 3 = inverted, 4 = left tilt)
    plxMat.Task.PercSgl       = nan(plxMat.NTrials,1);  % Percent trials w/ singleton distractor
    plxMat.Task.PercCatch     = nan(plxMat.NTrials,1);  % Percent catch trials

    if(size(plxMat.TrialMat_INFOcodes, 2) >= 34)
        plxMat.Task.DistFix = nan(plxMat.NTrials,1);  % fixed (1) or random (0) distractor color block
    end

    if(size(plxMat.TrialMat_INFOcodes, 2) >= 35)
        plxMat.Task.ProbCue = nan(plxMat.NTrials,1);  % is one location more often containing a target item?
    end
end

if(any(plxMat.TrialMat_EVcodes(:) == EV.Stimulation_))  % some form of stimulation ocurred (microstim, ultrasound)
    plxMat.Task.StimTrial = nan(plxMat.NTrials,1);
    plxMat.Task.StimEV    = nan(plxMat.NTrials,1);
    plxMat.Task.StimTime  = nan(plxMat.NTrials,1);
end

% loop over trials and identify relevant information
for(t=1:plxMat.NTrials)
    cInf = plxMat.TrialMat_INFOcodes(t,:);  % get the info codes corresponding to the current trial

    if(cInf(1) == 3000)
        plxMat.Task.TaskType{t} = 'MG';
    elseif(cInf(1) == 3999)
        plxMat.Task.TaskType{t} = 'Search';
    else
        plxMat.Task.TaskType{t} = 'NA';
    end

    % the numerical values that are subtracted are defined in the tempo
    % configuration files. Unfortunately this has to be hard coded and is
    % arbitrary in its definition. Be careful and double check thouroughy.

% position and angle assignment
%      ____________________
%     |  135 |  90  |   45 |
%     |  (7) |  (0) |  (1) |
%     |______|______|______|
%     |  180 |      |   0  |
%     |  (6) |  *   |  (2) |
%     |______|______|______|
%     |  225 |  270 |  315 |
%     |  (5) |  (4) |  (3) |
%     |______|______|______|
%
    switch char(plxMat.Task.TaskType(t))
        case 'MG'
            plxMat.Task.GoCue(t)        = plxMat.Task.FixSpotOff(t);
            plxMat.Task.TargetLoc(t)    = cInf(19) - 3000;
            plxMat.Task.MGcolor(t)      = cInf(20) - 3000;
            plxMat.Task.TargetCol_R(t)  = cInf(21) - 3000;
            plxMat.Task.TargetCol_G(t)  = cInf(22) - 3000;
            plxMat.Task.TargetCol_B(t)  = cInf(23) - 3000;

            plxMat.Task.Eccentricity(t) = cInf(24) - 3000;
            plxMat.Task.SOA(t)          = cInf(35) - 3000;
            plxMat.Task.TrialType(t)    = cInf(30) - 3000;

            plxMat.Task.Sess.FixSpotOffTEMPO(t) = plxMat.Task.StimOnset(t) + plxMat.Task.SOA(t);

            % PD event for fixation spot offset
            cpd = find(plxMat.TrialMat_PDtimes(t,:) >= plxMat.Task.FixSpotOffTEMPO(t), 1, 'first');
            if(~isempty(cpd))
                plxMat.Task.FixSpotOff(t) = plxMat.TrialMat_PDtimes(t,cpd);
            end

            plxMat.Task.MaxSaccDur(t)   = cInf(23) - 3000;
            plxMat.Task.MaxSaccTime(t)  = cInf(24) - 3000;
            plxMat.Task.TimeOut(t)      = cInf(25) - 3000;
            plxMat.Task.RewardDur(t)    = cInf(10) - 3000;
            plxMat.Task.RewardOffset(t) = cInf(11) - 3000;
            plxMat.Task.TargetHold(t)   = cInf(25) - 3000;
            plxMat.Task.Gains_EyeX(t)   = cInf(31) - 1000;
            plxMat.Task.Gains_EyeY(t)   = cInf(33) - 1000;
            plxMat.Task.Gains_XYF(t)    = cInf(32) - 1000;
            plxMat.Task.Gains_YXF(t)    = cInf(34) - 1000;

            plxMat.Task.TargetPos(t) = mod((360-plxMat.Task.TargetLoc(t)+90)/45,8); % transform angle to index position

        case 'Search'
            plxMat.Task.GoCue(t)         =  0;  % corresponds to search array 
            plxMat.Task.ArrayStruct(t)   =  cInf(2) - 4001;
            plxMat.Task.DistHomo(t)      =  cInf(3) - 4050;
            plxMat.Task.SingMode(t)      =  cInf(4) - 4060;
            plxMat.Task.SetSize(t)       =  cInf(5) - 4100;
            plxMat.Task.TargetType(t)    =  cInf(6) - 4150;
            plxMat.Task.TrialType(t)     =  cInf(7) - 4150;
            plxMat.Task.Eccentricity(t)  =  cInf(8) - 4250;
            plxMat.Task.Singleton(t)     = (cInf(9)- 1111)/1111;
            plxMat.Task.TargetLoc(t)     =  cInf(12) - 5000;
            plxMat.Task.DistLoc(t)       =  cInf(13) - 5500;
            plxMat.Task.IsCatch(t)       =  cInf(14) - 4300; % make sure code is 0 for no catch and 1 for catch
            plxMat.Task.SingleDistCol(t) =  cInf(15) - 4650;
            plxMat.Task.DistOri(t)       =  cInf(16) - 4660;
            plxMat.Task.TargetOri(t)     =  cInf(17) - 4670;
            plxMat.Task.PercSgl(t)       =  cInf(18) - 4700;
            plxMat.Task.PercCatch(t)     =  cInf(19) - 4800;
            plxMat.Task.BlockNum(t)      =  cInf(20) - 4900;
            plxMat.Task.SOA(t)           =  cInf(21) - 6000;

            plxMat.Task.MaxSaccDur(t)    = cInf(23) - 3000;
            plxMat.Task.MaxSaccTime(t)   = cInf(24) - 3000;
            plxMat.Task.TimeOut(t)       = cInf(25) - 3000;
            plxMat.Task.RewardDur(t)     = cInf(26) - 3000;
            plxMat.Task.RewardOffset(t)  = cInf(27) - 3000;
            plxMat.Task.TargetHold(t)    = cInf(28) - 3000;
            plxMat.Task.Gains_EyeX(t)    = cInf(30) - 1000;
            plxMat.Task.Gains_EyeY(t)    = cInf(31) - 1000;
            plxMat.Task.Gains_XYF(t)     = cInf(32) - 1000;
            plxMat.Task.Gains_YXF(t)     = cInf(33) - 1000;

            if(size(plxMat.TrialMat_INFOcodes, 2) >= 34)
                plxMat.Task.DistFix(t)   = abs(cInf(34) - 4682);
            end

            if(size(plxMat.TrialMat_INFOcodes, 2) >= 35)
                plxMat.Task.ProbCue(t)   = cInf(35) - 4690;
            end

            switch cInf(10)
                case 8888
                    plxMat.Task.TargetHemi(t) =  char('m');
                case 8200
                    plxMat.Task.TargetHemi(t) =  char('r');
                case 8100
                    plxMat.Task.TargetHemi(t) =  char('l');
            end

            switch cInf(11)
                case 8888
                    plxMat.Task.DistHemi(t) =  char('m');
                case 8200
                    plxMat.Task.DistHemi(t) =  char('r');
                case 8100
                    plxMat.Task.DistHemi(t) =  char('l');
            end

            if(plxMat.Task.SetSize(t) == 1) % with only one item it is a detection task
                plxMat.Task.TaskType{t} = 'Det';
            elseif(plxMat.Task.SingMode(t) == 1)  % double check this one
                plxMat.Task.TaskType{t} = 'Cap';
            elseif(size(plxMat.TrialMat_INFOcodes, 2) >= 35)
                if(cInf(35)-4690 > 0 && cInf(35)-4690 < 9)
                    plxMat.Task.TaskType{t} = 'ProbCue';
                end
            end

            plxMat.Task.TargetPos(t) = mod((360-plxMat.Task.TargetLoc(t)+90)/45,8); % transform angle to index position
            plxMat.Task.DistPos(t)   = mod((360-plxMat.Task.DistLoc(t)  +90)/45,8);
    end  % switch char(plxMat.Task.TaskType(t))

    % correct SRT time with the Go-signal event (PD event for fixation offset)
    plxMat.Task.SRT(t)     = plxMat.Task.Saccade(t) - plxMat.Task.GoCue(t);
%    plxMat.Task.SaccEnd(t) = plxMat.Task.SaccEnd(t) - plxMat.Task.StimOnsetToTrial(t);

    % check for the ocurrence of stimulation
    if(isfield(plxMat.Task, 'StimEV'))
        stimpos = find(plxMat.TrialMat_EVcodes(t,:) == 666);
        if(~isempty(stimpos))
            plxMat.Task.StimTrial(t) = 1;
            if(length(stimpos) > 1)
                warning(['More than one stimulation occurrence in trial ',int2str(t),'! Only log the first one and ignore the rest!']);
                stimpos(2:end) = [];
            end

            plxMat.Task.StimEV(t) = plxMat.TrialMat_EVtimes(t,stimpos);

            StimTime = plxMat.Task.StimEV(t) - plxMat.Task.StimOnsetToTrial(t);

            % reconstruct stimulation period
            switch char(plxMat.Task.TaskType(t))
                case 'MG'
                    if(StimTime < -50)
                        plxMat.Task.StimTime(t) = 1;
                    elseif(StimTime > -80 && StimTime < 0)
                        plxMat.Task.StimTime(t) = 2;
                    elseif(StimTime > 0 && StimTime < plxMat.Task.SOA(t) - 80)
                        plxMat.Task.StimTime(t) = 3;
                    elseif(StimTime > plxMat.Task.SOA(t) - 80 && StimTime < plxMat.Task.SOA(t))
                        plxMat.Task.StimTime(t) = 4;
                    elseif(plxMat.Task.Correct(t) == 0 && isnan(StimTime))
                        plxMat.Task.StimTrial(t) = 0;
                    else
                        warning('Something went wrong identifying stimulation period!');
                    end
                case 'Det'
                    if(StimTime < -50)
                        plxMat.Task.StimTime(t) = 1;
                    elseif(StimTime < 30 && StimTime > -30)
                        plxMat.Task.StimTime(t) = 2;
                    elseif(StimTime > 40 && StimTime < 100)
                        plxMat.Task.StimTime(t) = 3;
                    elseif(plxMat.Task.Correct(t) == 0 && isnan(StimTime))
                        plxMat.Task.StimTrial(t) = 0;
                    else
                        warning('Something went wrong identifying stimulation period!');
                    end
            end
        else
            % no stimulation in this trial
            plxMat.Task.StimTrial(t) = 0;
            plxMat.Task.StimTime(t)  = 0;
        end
    end
end  % for(t=1:plxMat.NTrials)

% Check for the presence of stimulation blocks in the intertrial interval

if(ITIstim == 1)
    % find code for stimulation end to determine blocks
    StimBlockEnd = find(plxMat.Event_Codes == EV.StimEnd_);

    if(~isempty(StimBlockEnd))  % data files contains flanking event codes for the stimulation block
        plxMat.Task.StimBlockEndTime = plxMat.Event_Times(StimBlockEnd);  % corresponding times

        ps = plxMat.Event_Codes == EV.StimStart_ | plxMat.Event_Codes == EV.Stimulation_;
        StimBlockstart = nan(1,length(StimBlockEnd));
        StimBlockCond  = nan(1,length(StimBlockEnd));
        plxMat.Task.StimTimes = [];

        for(s=1:length(StimBlockEnd))
            if(s==1)
                pps = plxMat.Event_Times <= plxMat.Task.StimBlockEndTime(s);
            else
                pps = plxMat.Event_Times <= plxMat.Task.StimBlockEndTime(s) & plxMat.Event_Times > plxMat.Task.StimBlockEndTime(s-1);
            end

            StimBlockstart(s) = find(ps == 1 & pps == 1, 1, 'first');

            if(plxMat.Event_Codes(StimBlockstart(s)) == EV.StimStart_)  % Non-Stim Block
                StimBlockCond(s) = 0;
            else % Stim Block
                StimBlockCond(s) = 1;
                stimtim = plxMat.Event_Times(ps == 1 & pps == 1 & plxMat.Event_Codes == EV.Stimulation_);
                plxMat.Task.StimTimes = [plxMat.Task.StimTimes; stimtim(2:end)];
            end
        end

        plxMat.Task.StimBlockStartTime = plxMat.Event_Times(StimBlockstart);
        plxMat.Task.StimBlockDur       = plxMat.Task.StimBlockEndTime - plxMat.Task.StimBlockStartTime;
        plxMat.Task.StimBlockPeriod    = median(diff(plxMat.Task.StimBlockStartTime));

        plxMat.Task.StimBlockTrial = nan(size(plxMat.Task.Correct));
        plxMat.Task.StimBlockCond  = nan(size(plxMat.Task.Correct));
        plxMat.Task.StimTrialTime  = plxMat.TrialStartEndDur(:,1);

        plxMat.Task.StimBlockCond  = nan(size(plxMat.Task.Correct));

        StimTrialTime  = plxMat.TrialStartEndDur(:,1);
        for(s=1:length(StimBlockstart))
            if(s<length(StimBlockstart))
                cspos = plxMat.TrialStartEndDur(:,1) > plxMat.Task.StimBlockEndTime(s) & plxMat.TrialStartEndDur(:,1) < plxMat.Task.StimBlockStartTime(s+1);
            else
                cspos = plxMat.TrialStartEndDur(:,1) > plxMat.Task.StimBlockEndTime(s);
            end
            StimTrialTime(cspos)  = StimTrialTime(cspos) - plxMat.Task.StimBlockEndTime(s);
            plxMat.Task.StimBlockTrial(cspos) = s;
            plxMat.Task.StimBlockCond(cspos)  = StimBlockCond(s);
            plxMat.Task.StimTrialTime(cspos)  = plxMat.Task.StimTrialTime(cspos) - plxMat.Task.StimBlockEndTime(s); % plxMat.Task.StimTimes(LongGapPos(s)+1);
        end

    else % no flanking codes in the data file, try to reconstruct.
        plxMat.Task.StimTimes = plxMat.Event_Times(plxMat.Event_Codes == EV.Stimulation_); % find event codes for stimulation
        StimGap = diff(plxMat.Task.StimTimes);  % determine spacing between subsequent stimulations
        StimBlockstart = [0; find(StimGap-median(plxMat.Task.StimGap) > mad(StimGap))]+1; % Blocks are identified by a long interval preceding the stimulation

        plxMat.Task.StimBlockStartTime = plxMat.Task.StimTimes(StimBlockstart);  % corresponding times

        plxMat.Task.StimBlockDur = nan(length(StimBlockstart),1); % determine duration for the stimulation block

        for(s=1:length(StimBlockstart)-1)
            plxMat.Task.StimBlockDur(s) = diff(plxMat.Task.StimTimes([StimBlockstart(s),StimBlockstart(s+1)-1]));
        end
        plxMat.Task.StimBlockDur(s+1) = plxMat.Task.StimTimes(StimBlockstart(end)) - plxMat.Task.StimTimes(StimBlockstart(s));

        plxMat.Task.StimBlockPeriod = median(plxMat.Task.StimBlockDur);  % period of two subsequent stimulation blocks (approximate)

        StimTrialTime  = plxMat.TrialStartEndDur(:,1);
        for(s=1:length(StimBlockstart)-1)
            cspos = plxMat.TrialStartEndDur(:,1) > plxMat.Task.StimBlockStartTime(s) & plxMat.TrialStartEndDur(:,1) < plxMat.Task.StimBlockStartTime(s+1);
            StimTrialTime(cspos)  = StimTrialTime(cspos) - plxMat.Task.StimBlockStartTime(s);
        end
        cspos = plxMat.TrialStartEndDur(:,1) > plxMat.Task.StimBlockStartTime(s+1);
        StimTrialTime(cspos)  = StimTrialTime(cspos) - plxMat.Task.StimBlockStartTime(s+1);

        % identify long ITI corresponding to non-stimulation pauses
        LongGapPos = find(diff(plxMat.TrialStartEndDur(:,1)) > 0.8*median(plxMat.Task.StimBlockDur));

        plxMat.Task.StimBlockTrial = nan(size(plxMat.Task.Correct));
        plxMat.Task.StimBlockCond  = nan(size(plxMat.Task.Correct));
        plxMat.Task.StimTrialTime  = plxMat.TrialStartEndDur(:,1);

        for(s=1:length(LongGapPos))
            if(s<length(LongGapPos))
                cspos = LongGapPos(s)+1 : LongGapPos(s+1);
            else
                cspos = LongGapPos(s)+1 : plxMat.NTrials;
            end

            % tag each trial with the number of the preceding stimulation block and determine time relative to stimulation
            plxMat.Task.StimBlockTrial(cspos) = s;
            plxMat.Task.StimTrialTime(cspos)  = plxMat.Task.StimTrialTime(cspos) - plxMat.TrialStartEndDur(LongGapPos(s)+1,1); % plxMat.Task.StimTimes(LongGapPos(s)+1);

            if(StimTrialTime(LongGapPos(s)) > StimTrialTime(LongGapPos(s)+1))
                plxMat.Task.StimBlockCond(cspos) = 0;
            else
                plxMat.Task.StimBlockCond(cspos) = 1;
            end
        end
    end

    if(keep_sess == 0) % clean up some data that is not wanted
        plxMat = rmfield(plxMat, 'Event_Codes');
        plxMat = rmfield(plxMat, 'Event_Times');
    end
end

% ad hoc hac to determine location with higher probability to have a target

LocCnt_all = hist(plxMat.Task.TargetPos,-0.5:7);
plxMat.Task.LocProb_all = LocCnt_all/sum(LocCnt_all);

LocCnt_hit = hist(plxMat.Task.TargetPos(plxMat.Task.Correct),-0.5:7);
plxMat.Task.LocProb_hit = LocCnt_hit/sum(LocCnt_hit);

plxMat.Task.LocErr = 1 - (LocCnt_hit./LocCnt_all);

plxMat.Task.CuePos = find(LocCnt_hit > 2*mean(LocCnt_hit),1,'first');
if(isempty(plxMat.Task.CuePos))
    plxMat.Task.CuePos = NaN;
end

% ____________________________________________________________________________%
%% define error outcome
% This needs to be checked and coded more dynamically
plxMat.Task.error_names = {'False', 'Early', 'Late', 'FixBreak', 'HoldError', ...
    'CatchErrorGo', 'CatchErrorNoGo'};

plxMat.Task.error = nan(plxMat.NTrials,1);
plxMat.Task.error(plxMat.Task.Correct == 1) = 0;

false_resp = logical(sum(bsxfun(@eq, plxMat.TrialMat_EVcodes, EV.Error_sacc),2));
if(any(false_resp))
    plxMat.Task.error(false_resp) = 1;
end

early_resp = logical(sum(bsxfun(@eq, plxMat.TrialMat_EVcodes, EV.EarlySaccade_),2));
if(any(early_resp))
    plxMat.Task.error(early_resp) = 2;
end

fix_break = logical(sum(bsxfun(@eq, plxMat.TrialMat_EVcodes, EV.FixError_),2));
if(any(fix_break))
    plxMat.Task.error(fix_break) = 4;
end

hold_err = logical(sum(bsxfun(@eq, plxMat.TrialMat_EVcodes, EV.BreakTFix_),2));
if(any(hold_err))
    plxMat.Task.error(hold_err) = 5;
end

catch_go = logical(sum(bsxfun(@eq, plxMat.TrialMat_EVcodes, EV.CatchIncorrectG_),2));
if(any(catch_go))
    plxMat.Task.error(catch_go) = 6;
end

catch_hold = logical(sum(bsxfun(@eq, plxMat.TrialMat_EVcodes, EV.CatchIncorrectG_),2));
if(any(catch_hold))
    plxMat.Task.error(catch_hold) = 7;
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
        elseif(any(plxMat.TrialMat_EVcodes(cet,:) == EV.Abort_))
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
%% get EEG channels
if(any(ADchan.is_eeg))
    plxMat.EEG_names = {};
    plxMat.EEGchan   = [];
    eeg_probes = fieldnames(ADchan.EEG);

    for(e = 1:length(eeg_probes))
        cprobe = cell2mat(eeg_probes(e));
        eeg_chan = ADchan.EEG.(cprobe);

        cAD = sprintf('AD%02d', eeg_chan);
%        if(~isfield(plxMat, cAD) && ~isfield(plxMat.EEG, cAD))
        if(~isfield(plxMat, cAD))
            warning(['Did not find a corresponding analog channel for the EEG probe at ',cprobe,'']);
        else
            % pre-allocate
            plxMat.EEG.(cprobe) = AD_mat;

            % just another painful trial loop...
            for(ct = 1:plxMat.NTrials)
                cpos = find(timebin_mat(ct,:) >= PreOnsetTime);
                [~, spos] = intersect(plxMat.timevec, timebin_mat(ct,cpos));

                if(~isempty(cpos)) % discard data from trials without stimulus onset
                    plxMat.EEG.(cprobe)(ct,spos) = plxMat.(cAD)(ct,cpos);
                end
            end
            plxMat = rmfield(plxMat, cAD);
            %             plxMat.EEG.([cprobe,'_z']) = plxMat.([cAD,'_z']);
            %             plxMat   =  rmfield(plxMat, [cAD,'_z']);

            chanpos = find(plxMat.ADChannel == eeg_chan - 1);
            plxMat.EEG.([cprobe,'_info']).orig      = plxMat.ADName(    chanpos);
            plxMat.EEG.([cprobe,'_info']).gain      = plxMat.ADGain(    chanpos);
            plxMat.EEG.([cprobe,'_info']).pregain   = plxMat.ADPreGain( chanpos);
            plxMat.EEG.([cprobe,'_info']).channel   = plxMat.ADChannel( chanpos);
            plxMat.EEG.([cprobe,'_info']).rate      = plxMat.ADSmplRate(chanpos);
            plxMat.EEG.([cprobe,'_info']).timestamp = plxMat.ADStartRec(chanpos); % this should be corrected by now

            plxMat.EEG_names = [plxMat.EEG_names, cprobe];
            plxMat.EEGchan   = [plxMat.EEGchan, eeg_chan];

            if(isfield(plxMat,[cAD,'_session']))
                if(keep_sess == 1)
                    plxMat.EEG.([cprobe,'_session']) = plxMat.([cAD,'_session']);
                end
                plxMat = rmfield(plxMat, [cAD,'_session']);
            end
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

%% get channel for stimulation
% stimulation time course
if(~isempty(ADchan.Stim1) || ~isempty(ADchan.Stim2))
    if(isfield(plxMat, sprintf('AD%02d', ADchan.Stim1)))
        eS = sprintf('AD%02d', ADchan.Stim1);
        stim_chan = ADchan.Stim1;
    elseif(isfield(plxMat, sprintf('AD%02d', ADchan.Stim2)))
        eS = sprintf('AD%02d', ADchan.Stim2);
        stim_chan = ADchan.Stim2;
    else
        stim_chan = [];
    end
else
    stim_chan = [];
end

if(~isempty(stim_chan))
    % pre-allocate
    plxMat.Stim.Pulse = AD_mat;

    % just another painful trial loop...
    for(ct = 1:plxMat.NTrials)
        cpos = find(timebin_mat(ct,:) >= PreOnsetTime);
        [~, spos] = intersect(plxMat.timevec, timebin_mat(ct,cpos));

        if(~isempty(cpos)) % discard data from trials without stimulus onset
            plxMat.Stim.Pulse(ct,spos) = plxMat.(eS)(ct,cpos);
        end
    end
    plxMat = rmfield(plxMat, eS);
    %         plxMat =  rmfield(plxMat, [eX,'_z']);

    chanpos = find(plxMat.ADChannel == stim_chan - 1);
    plxMat.Stim.info.orig      = plxMat.ADName(    chanpos);
    plxMat.Stim.info.gain      = plxMat.ADGain(    chanpos);
    plxMat.Stim.info.pregain   = plxMat.ADPreGain( chanpos);
    plxMat.Stim.info.channel   = plxMat.ADChannel( chanpos);
    plxMat.Stim.info.rate      = plxMat.ADSmplRate(chanpos);
    plxMat.Stim.info.timestamp = plxMat.ADStartRec(chanpos); % this should be corrected by now
end

% ____________________________________________________________________________%
%% Check if there are remaining AD channels that were not interpreted
% there must be a better way to find fieldnames according to a pattern and check
% if all analog channels were renamed. For now this seems to be at least a
% working solution.
leftfields = char(fieldnames(plxMat));
lefty      = sort(leftfields,2);
leftovers  = strcmp('AD',lefty(:,end-1:end));
if(sum(leftovers) > 0)
    warning('Not all analog channels could be interpreted. Check carefully!')
    disp('leftover analog channels:');
    disp(leftfields(leftovers,:));
else
    % clean up
    plxMat = rmfield(plxMat, 'ActAD');
    plxMat = rmfield(plxMat, 'ADName');
    plxMat = rmfield(plxMat, 'ADGain');
    plxMat = rmfield(plxMat, 'ADPreGain');
    plxMat = rmfield(plxMat, 'ADChannel');
    plxMat = rmfield(plxMat, 'ADSmplRate');
    plxMat = rmfield(plxMat, 'ADStartRec');
end

% ____________________________________________________________________________%
%% Get more detailed description of eye traces
plxMat = PLX_EyeTrace(plxMat);

% ____________________________________________________________________________%
%% Determine data size
if(exist('ByetSize','file'))
    fprintf('Resulting structure has a size of %s', ByetSize(plxMat));
end

% ____________________________________________________________________________%
%% helper functions

function ctm = PLXin_get_event_time(EVcode, ctrial)

% use of evalin usually not recommended. Here, as kind of inline code
EVc = evalin('caller',['plxMat.TrialMat_EVcodes(',int2str(ctrial),',:)']);
EVt = evalin('caller',['plxMat.TrialMat_EVtimes(',int2str(ctrial),',:)']);
ctm = NaN;

p = EVc == EVcode;
if(sum(p)==1)
    ctm = EVt(p);
elseif(sum(p)>1)
    warning('Something weird occurred: more than one event found. Check carefully!');
end



