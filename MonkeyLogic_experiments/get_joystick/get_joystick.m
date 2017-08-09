% get_joystick (timing script)
%
% Initial timing file for teaching the NHP to use a joystick in order to get a reward.
%
% A stimulus is shown and the monkey has to pull the joystick within a specified period of <wait_start>.
% If <PullRew> is set to 1, a reward will be given for pulling the lever. After a variable time interval
% in the range of <dimmMin> and <dimmMax> the stimulus is changing (as defined in the condition file) and
% the monkey has to release the joystick to receive another reward (incremental scheme according to number
% of hits or for consecutive correct trials).
%
% Be aware: the current implementation expects a reactive behaviour to the onset of the target object
% to start a trial by pressing the joystick, hold it for the course of the trial and release it as
% a stimulus change. This will later on result in the requirement to fixate and pull simultaneous
% There are alternative approaches that might be more desired, for example by waiting for a joystick
%  and show the target object as result of the press, and again wait for a release as task related
% first the press is required later on, followed by fixation.
% Another option might be to first show the target object, wait for a fixation to start the trial,
% a brief lever press as task related response (press+release within a short time window).
%
% Also be aware, that this timing file is made to be used with a analogue joystick as response device.
% The use of an analogue joystick changes the calls and logic of the eyejoytrack function.
%
% wolf zinke, Feb. 2014

%% ### ToDo: ### %%
    % use the gen function in the condition file to create the target item and make its features (size, luminance) editable.
    % define ITI as variable in TrialRecord (done!) and use the ITI end time instead of ML iti mechanism.
    % get subject name and paradigm from MLConfig (not accessible from timing file right now).

eventmarker(15);  %  start pre-trial

scriptversion='20140425'   % update this information to control if MonkeyLogic is really using an updated timing file

%% Define Trial Variables
tblpath = pwd;  % path to directory for trial table files
cNHP   = 'NHP_name';        % MLConfig.SubjectName    - MLConfig could not be accessed from the timing file
cPRDGM = 'get_joystick';  % MLConfig.ExperimentName - MLConfig could not be accessed from the timing file

% set editable variables
editable('Jpos0x', 'Jpos0y', 'JoyPullRad', 'JoyRelRad', 'dimmMin', 'dimmMax', 'dimmStep', 'wait_resp', 'max_resp', 'time_out', 'pull_pause', 'minITI', 'maxITI', 'ITIstep','PullRew', 'RewInc2P', 'RewInc3P', 'RewInc4P', 'RewInc5P');

% initialize the random number generator
rng('shuffle', 'twister');

% reward times
PullRew    =     0;    % Give reward for pulling the joystick
rew_dur    =    50;    % duration of reward pulse (not to be changed!)
rew_gap    =   200;    % gap between two reward pulses (period = rew_dur+rew_gap) (not to be changed!)
rew_lag    =   100;    % introduce a brief delay prior reward
RewIncConsecutive = 1; % increase reward for subsequent correct trials. Otherwise reward will increase with the number of hits

% set center position of joystick
Jpos0x     =   0;      % zero position X
Jpos0y     =   0;      % zero position Y
JoyPullRad = 1.5;      % threshold to detect a elevation of the joystick
JoyRelRad  =   1;      % threshold to detect the release of the joystick

% dimming parameters
dimmMin    =  500;     % minimum time to dimming
dimmMax    = 2000;     % maximum time to dimming
dimmStep   =   10;     % steps of possible dimming times

% inter trial times (replace this with a more sophisticated function)
minITI  =   500;      % minimum range of the time period that separates subsequent presentation of the stimulus/trials
maxITI  =  3000;      % maximum range of the time period that separates subsequent presentation of the stimulus/trials
ITIstep =   100;
ITIvec  = minITI : ITIstep : maxITI;     % possible ITI's
cITI    = ITIvec(randi(length(ITIvec))); % ITI used for the current trial

% joystick response parameters
wait_start = 2000;    % maximal time period to press the lever in order to start a trial.
wait_resp  =  100;    % valid response only accepted after this initial period
max_resp   = 1000;    % maximal time accepted for a valid responses

time_out   = 1000;    % set wait period for bad monkeys (do not combine with jittered ITIs)
pull_pause = 6000;    % this is the minimum time passed before a trial starts after random lever presses

%% Assign stimulus items
% assign items specified in the condition file to meaningful names
% numbers correspond to the number following TaskObject# in the condition file
JoyZero     =  1;     % this might be a simple trick to lock the joystick to a zero position by using a Fixation object
StimInitial =  2;     % first item state
StimChange  =  3;     % item after contrast change

%% initialize trial variables
on_track   =   0;     % use this flag to indicate occurrences of errors in the trial
rt         = NaN;
RelTime    = NaN;
RewTime    = NaN;
WaitStart  = NaN;
TrialZero  = NaN;
TstartEff  = NaN;
TRelEff    = NaN;
DimmOnTime = NaN;
DimmEff    = NaN;
FixOn      = NaN;
FixOff     = NaN;

TrialRecord.Tstart(TrialRecord.CurrentTrialNumber) = NaN;
TrialRecord.Tend(  TrialRecord.CurrentTrialNumber) = NaN;

%% determine number of reward pulses (incremental reward scheme)
if(~isfield(TrialRecord, 'CorrCount'))
    TrialRecord.CorrCount = 0;
end

if(RewIncConsecutive == 1)
   % increase reward for consecutive correct trials
    RewInc2P = 1;   %   two pulses
    RewInc3P = 2;   % three pulses
    RewInc4P = 3;   %  four pulses
    RewInc5P = 4;   %  five pulses (jackpot)

    if(TrialRecord.CorrCount < RewInc2P-1)
        rew_Npulse = 1;
    elseif(TrialRecord.CorrCount < RewInc3P-1)
        rew_Npulse = 2;
    elseif(TrialRecord.CorrCount < RewInc4P-1)
        rew_Npulse = 3;
    elseif(TrialRecord.CorrCount < RewInc5P-1)
        rew_Npulse = 4;
    else
        rew_Npulse = 5;
    end
else
    cNumHit = sum(TrialRecord.TrialErrors == 0);

    % increase rewards after a defined number of trials was achieved
    RewInc2P =  75;   %   two pulses
    RewInc3P = 200;   % three pulses
    RewInc4P = 350;   %  four pulses
    RewInc5P = 450;   %  five pulses (jackpot)

    if(cNumHit < RewInc2P)
        rew_Npulse = 1;
    elseif(cNumHit < RewInc3P)
        rew_Npulse = 2;
    elseif(cNumHit < RewInc4P)
        rew_Npulse = 3;
    elseif(cNumHit < RewInc5P)
        rew_Npulse = 4;
    else
        rew_Npulse = 5;
    end
end

% split the stimulus change/dimming times across conditions to allow for some control
num_cnd = length(unique(TrialRecord.ConditionsThisBlock)); % This variable is behaving strange and not as I understand its role. The use of unique ensures that it shows the actual number of conditions used in this block!
ccnd    = TrialRecord.CurrentCondition

itvtm   = linspace(dimmMin, dimmMax, num_cnd+1);
Dimmvec = dimmMin : dimmStep : dimmMax; % possible stimulus change times
Dimmvec = Dimmvec(Dimmvec >= itvtm(ccnd) & Dimmvec < itvtm(ccnd+1)); % possible stimulus change times for current condition

cDimm   = Dimmvec(randi(length(Dimmvec)))        % stimulus change time used for the current trial

%%%%###########################################%%%%
%%%%########  Ensure joystick release   #######%%%%
%%%%###########################################%%%%
%% check for joystick position.
% If the joystick is pulled, a defined period will be waited for a release,
% and this trial is aborted either after this period passed or the lever is released.
% This trial restart is necessary to avoid ML crashes due too long trial durations.

reposition_object(JoyZero, Jpos0x, Jpos0y);  % correct for offsets of the joystick analogue output

if(~isfield(TrialRecord, 'Jreleased'))
    TrialRecord.Jreleased = 1;   % track the joystick state across trials
end

% check if lever remains released
[JoyRel PullTime] = eyejoytrack('holdtarget', JoyZero, JoyRelRad, 100);

if(~JoyRel | TrialRecord.Jreleased == 0)
%      if(cDist > JoyRelRad)
    if(~JoyRel) % wait for release
        TrialRecord.Jreleased = 0;
        [JoyRel] = eyejoytrack('acquiretarget', JoyZero, JoyRelRad, pull_pause);

        if(JoyRel)
            eventmarker(39); % lever released
        end

    else % make sure that joystick is not pulled for a sufficient long time before starting a new trial
        [JoyRel PullTime] = eyejoytrack('holdtarget', JoyZero, JoyRelRad, pull_pause);

        if(JoyRel)
            TrialRecord.Jreleased = 1;
        end
    end

    TrialRecord.CorrCount = 0;
    trialerror(8);
    set_iti(0);
else
    on_track = 1;
    TrialRecord.Jreleased = 1;
end

eventmarker(16); % End pre trial

%%%%###########################################%%%%
%%%%##########  Start the Trial now   #########%%%%
%%%%###########################################%%%%

if(on_track)
%% #### start Trials #### %%
    % show fix spot to encourage starting the trial
    eventmarker(10);  % end ITI

    if(~isfield(TrialRecord, 'NextTrial'))
        TrialRecord.NextTrial = -1;
    end

    StimOnTime = toggleobject(StimInitial, 'EventMarker', 35, 'Status', 'on');

    % get the computer time for the trial start as well
    TrialRecord.Tstart(TrialRecord.CurrentTrialNumber) = toc(uint64(1));

    if(TrialRecord.CurrentTrialNumber > 1)
        ITIeff = TrialRecord.Tstart(TrialRecord.CurrentTrialNumber) - TrialRecord.Tend(TrialRecord.CurrentTrialNumber-1);
    else
        ITIeff = NaN;
    end

%% #### Wait Trial Start #### %%
    % wait for joystick touch
    % this seems counter-intuitive first, but 'holdtarget' is interpreted here as being in zero position.

    eventmarker(11);  % wait for start
    [JoyRel PullTime] = eyejoytrack('holdtarget', JoyZero, JoyPullRad, wait_start);
    eventmarker(12);  % end waiting for start
    evobj.CodeTimes = trialtime;

    if(JoyRel && on_track)
        ErrCode  = 9; % not started
        StimOffTime = toggleobject(StimInitial, 'EventMarker', 36,'Status', 'off');
        on_track = 0;
    elseif(on_track)
        %eventmarker([1 38]);  % lever pressed
        eventmarker(1);   % trial started
        eventmarker(38);  % lever pressed
        TrialRecord.Jreleased = 0;

        % get the computer time for the trial start as well
        TstartEff = toc(uint64(1));

        TrialZero = evobj.CodeTimes(1);
        WaitStart = TrialZero - StimOnTime;
        FixOn     = StimOnTime - TrialZero;
        DimmTime  = TrialZero  + cDimm;

        if(PullRew == 1)
            eventmarker(19);  % start pause
            idle(rew_lag);
            %eventmarker([20 96]);
            eventmarker(20);   % end pause
            eventmarker(96);   % reward start
            goodmonkey(rew_dur, 'NumReward', 1, 'PauseTime', rew_gap);
            eventmarker(96);   % reward end
       end

        [JoyRel ReleaseTime] = eyejoytrack('acquiretarget', JoyZero, JoyRelRad, DimmTime-trialtime);  % wait before dimming
    end

%% #### Wait Dimming #### %%
    % do not accept early responses
    if(JoyRel && on_track)
        eventmarker(39);  % lever released

        evobj.CodeTimes = trialtime;
        % get the computer time for the lever release as well
        TRelEff = toc(uint64(1));
        RelTime =  evobj.CodeTimes(1) - TrialZero;

        TrialRecord.Jreleased = 1;
        ErrCode  = 5;
        on_track = 0;
        StimOffTime = toggleobject(StimInitial , 'EventMarker', 36,'Status', 'off');
    elseif(on_track)
        DimmOnTime = toggleobject([StimInitial StimChange], 'EventMarker', 37);  % dimm target
        DimmOnTime = DimmOnTime - TrialZero;

        DimmEff = toc(uint64(1));
    end

%% #### Wait Response #### %%
    % wait for release of joystick
    if(on_track)
        [JoyRel ReleaseTime] = eyejoytrack('acquiretarget', JoyZero, JoyRelRad, max_resp);
        % get the computer time for the lever release as well

        if(~JoyRel)
            ErrCode  = 1; % no response
            on_track = 0;
            StimOffTime = toggleobject(StimChange, 'EventMarker', 36,'Status', 'off');
        else
            eventmarker(39);  % lever released

            evobj.CodeTimes = trialtime;
            RelTime = evobj.CodeTimes(1) - TrialZero;

            TRelEff = toc(uint64(1));

            TrialRecord.Jreleased = 1;

            if(ReleaseTime < wait_resp)
                ErrCode  = 5; % early response
                on_track = 0;
                StimOffTime = toggleobject(StimChange , 'EventMarker', 36,'Status', 'off');
            end
        end
    end

%% #### Trial End #### %%
    %eventmarker([2 17]);  % end of trial
    eventmarker(2);  % end of trial
    eventmarker(17);  % start post trial

    % finalize Trial
    if(on_track)  % correct trials
        ErrCode = 0;
        trialerror(0); % correct
        rt = ReleaseTime;  % double check this one (use text table and compare to RelTime-DimmOnTime!)
        TrialRecord.CorrCount = TrialRecord.CorrCount + 1;

        eventmarker(19);   % start pause
        idle(rew_lag);
        %eventmarker([96 20]);
        eventmarker(20);   % end pause

        StimOffTime = toggleobject(StimChange, 'EventMarker', 36,'Status', 'off');
        eventmarker(96);    % reward start
        evobj.CodeTimes = trialtime;

        goodmonkey(rew_dur, 'NumReward', rew_Npulse, 'PauseTime', rew_gap);
        eventmarker(96); % use this event twice to get an estimate of the reward duration, i.e. number of pulses.
        RewTime = evobj.CodeTimes(1) - TrialZero;
    else   % error occurred
        trialerror(ErrCode);
        rew_Npulse = NaN;
        TrialRecord.CorrCount = 0;
        cITI = cITI + time_out;   % add a punish time to the current ITI
    end

    FixOff = StimOffTime - TrialZero;

    % get the computer time for the trial end as well
    TrialRecord.Tend(TrialRecord.CurrentTrialNumber) = toc(uint64(1));

%%%%###########################################%%%%
%%%%##########   Stop the Trial now   #########%%%%
%%%%###########################################%%%%
    % set ITI
    %eventmarker([18 9]); % start ITI
    eventmarker(18); % end post trial
    eventmarker(9);  % start ITI
    set_iti(cITI);
    TrialRecord.NextTrial = TrialRecord.Tend(TrialRecord.CurrentTrialNumber) + cITI/1000;

    if(ErrCode ~= 8 && ErrCode ~= 9)
    %% write trial table
    % keep redundant information to check data and have some information as back up to reconstruct trials
    % Further, this text file should give fast and easy access to the behavioral performance. It is in a
    % format that allows to be read into R with the read.table command (use 'header=TRUE' option).
        cdt   = datestr(now,'yyyy_mm_dd');
        tblnm = fullfile(tblpath, [cNHP,'_', cPRDGM,'_', cdt, '.dat']);

        if(exist(tblnm) ~= 2 )  % create the table file
            tblptr = fopen(tblnm, 'w');
            fprintf(tblptr,'Date  Subject  Experiment  TFixOn  Tstart  TRel  Tend  ITIeff  TrialNo  BlockNo  CondNo  Result  WaitStart  ITI  NumRew  RT  TrialZero  FixOn  FixOff  DimmEff  DimmOn  DimmIntend  ReleaseTime  RewTime\n');
        else
            tblptr = fopen(tblnm, 'a');
        end

        fprintf(tblptr,'%s  %s  %s  %.4f  %.4f  %.4f  %.4f  %6d  %6d  %4d  %4d  %4d  %8d  %8d  %8d  %8d  %8d  %8d  %8d  %.4f  %8d  %8d  %6d  %6d \n',...
                cdt, cNHP, cPRDGM, TrialRecord.Tstart(TrialRecord.CurrentTrialNumber), TstartEff, TRelEff, ...
                TrialRecord.Tend(TrialRecord.CurrentTrialNumber), ITIeff, TrialRecord.CurrentTrialNumber,  ...
                TrialRecord.CurrentBlock, ccnd, ErrCode, WaitStart, cITI, rew_Npulse, ...
                rt, TrialZero, FixOn, FixOff, DimmEff, DimmOnTime, cDimm, RelTime, RewTime );
        fclose(tblptr);
    end
end

