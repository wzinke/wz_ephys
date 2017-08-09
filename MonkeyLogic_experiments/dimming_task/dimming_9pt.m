% dimming (timing script)
%
% Initial timing file for teaching the NHP to use a joystick in order to get a reward.
%
% This timing file is prepared to run with the GitHub version of
% MonkeyLogic (https://github.com/Dfreedmanlab/MonkeyLogic_stable). It is possible that it will cause
% errors when used with previous versions of MonkeyLogic.
%
% This task utilizes two Stimuli. One of them signals the animal not to pull the joystick (NoGo Stimulus),
% whereas the other stimulus (Go Stimulus) requires a pull.
%
% The NoGo stimulus is shown for a random interval in the range of <dimmMin> and <dimmMax>.  The animal has to
% pull the joystick within a time window between <wait_resp> and <max_resp> following the onset of the Go stimulus
% to receive a reward.
%
% Be aware, that this timing file is made to be used with a analogue joystick as response device.
% The use of an analogue joystick changes the calls and logic of the eyejoytrack function.
%
% wolf zinke, May. 2014

%% ### ToDo: ### %%
    % use the gen function in the condition file to create the target item and make its features (size, luminance) editable.
    % define ITI as variable in TrialRecord (done!) and use the ITI end time instead of ML iti mechanism.
    % get subject name and paradigm from MLConfig (not accessible from timing file right now).

eventmarker(15);  %  start pre-trial

%% Define Trial Variables
cdt   = datestr(now,'yyyy_mm_dd');
tblpath = ['C:/Users/user/Desktop/Training_Data/I_34/', cdt];  % path to directory for trial table files
cNHP   = 'I_34';      % MLConfig.SubjectName    - MLConfig could not be accessed from the timing file
cPRDGM = 'dimming';   % MLConfig.ExperimentName - MLConfig could not be accessed from the timing file

% set editable variables
editable('Jpos0x', 'Jpos0y', 'JoyPullRad','FixGrid', 'ForceLoc', 'LocRep', 'AccAll', 'FixStep', 'sqrsz', 'NGsqrcol', 'Gsqrcol', 'dimmMin', 'dimmMax', 'dimmStep', 'wait_resp', 'max_resp', 'wait_late', 'wait_rel', 'time_out', 'pull_pause', 'minITI', 'maxITI', 'ITIstep', 'BlockDur', 'BlockBreakMin', 'BlockBreakMax', 'useAddITI', 'RewInc2P', 'RewInc3P', 'RewInc4P', 'RewInc5P');

% initialize the random number generator
rng('shuffle', 'twister');

% reward times
rew_dur    =    50;    % duration of reward pulse (not to be changed!)
rew_gap    =   200;    % gap between two reward pulses (period = rew_dur+rew_gap) (not to be changed!)
rew_lag    =   100;    % introduce a brief delay prior reward
RewIncCons =     1;    % increase reward for subsequent correct trials. Otherwise reward will increase with the number of hits

% set centre position of joystick
Jpos0x     = -0.5;     % zero position X
Jpos0y     =  0.5;     % zero position Y
JoyPullRad =  0.5;    % threshold to detect a elevation of the joystick

% dimming parameters
dimmMin    = 1000;     % minimum time to dimming
dimmMax    = 2500;     % maximum time to dimming
dimmStep   =   10;     % steps of possible dimming times

NGsqrcol   = 0.80;     % grey value for pre-dimm item
Gsqrcol    = 0.50;     % grey value for post-dimm item
sqrsz      =    1;     % stimulus size (square)

TrialRecord.sqrszsize = sqrsz;
TrialRecord.NGsqrcol  = NGsqrcol;
TrialRecord.Gsqrcol   = Gsqrcol;

% stimulus position parameters
FixGrid    =    9;     % number of possible locations for the dimming item. <1> means only presented in centre
                       % other options are <5> and <9>
LocRep     =    2;     % how many trials are shown per location in a Block
AccAll     =    0;     % accept all results, i.e. hit/early/late to complete current trial.
                       % If 0, only correct trials will be accepted
FixStep    =    5;     % y and x distance of possible item locations
ForceLoc   =    0;     % if this number is between 1:9 it forces the stimulus location to that position


% inter trial times (replace this with a more sophisticated function)
minITI     =  1000;    % minimum time period between two subsequent stimulus presentation
maxITI     =  3500;    % maximum time period between two subsequent stimulus presentation
ITIstep    =    50;

ITIvec  = minITI : ITIstep : maxITI;     % possible ITI's
cITI    = ITIvec(randi(length(ITIvec))); % ITI used after the current trial (uniform distribution)

useAddITI = 0;

% define ITI according to a shifted and truncated exponential distribution
%  mu = (minITI + maxITI)/2
%  max_range = exp(-(minITI/mu));
%  min_range = exp(-(maxITI/mu));
%  cITI = ITIstep * round( (-mu*reallog((max_range-min_range)*rand(1,1)+min_range)) / ITIstep);

BlockDur      =  12;   % Duration of a training block in minutes
BlockBreakMin = 1.5;   % minimum break between two blocks in minutes
BlockBreakMax =   3;   % maximum break between two blocks in minutes

% joystick response parameters
wait_resp  =  150;     % valid response only accepted after this initial period
max_resp   = 1500;     % maximal time accepted for a valid responses

time_out   = 2500;     % set wait period for bad monkeys (do not combine with jittered ITIs)
pull_pause = 5000;     % this is the minimum time passed before a trial starts after random lever presses
wait_rel   = 1000;     % wait for release of joystick after response
wait_late  = 1000;     % just wait to check if a late response occurs.
%% Assign stimulus items
% assign items specified in the condition file to meaningful names
% numbers correspond to the number following TaskObject# in the condition file
JoyZero  =  1;         % this might be a simple trick to lock the joystick to a zero position by using a Fixation object
StimGo   =  2;         % stimulus that requires a response
StimNoGo =  3;         % stimulus that indicates no response

%% initialize trial variables
on_track      =   0;  % use this flag to indicate occurrences of errors in the trial
rt            = NaN;
TrialZero     = NaN;
TRelEff       = NaN;
FixOn         = NaN;
FixOff        = NaN;
FixOffEff     = NaN;
PullTime      = NaN;
TPullEff      = NaN;
RelTime       = NaN;
StimOffTime   = NaN;
NoGoOn        = NaN;
GoOn          = NaN;
NoGoEff       = NaN;
GoEff         = NaN;
RewEff        = NaN;
ITIeff        = NaN;

TrialRecord.Tstart(TrialRecord.CurrentTrialNumber) = NaN;
TrialRecord.Tend(  TrialRecord.CurrentTrialNumber) = NaN;

%% determine number of reward pulses (incremental reward scheme)
if(~isfield(TrialRecord, 'CorrCount'))
    TrialRecord.CorrCount = 0;
end

if(RewIncCons == 0)
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
else
   % increase reward for consecutive correct trials
    RewInc2P = 2;   %   two pulses
    RewInc3P = 3;   % three pulses
    RewInc4P = 4;   %  four pulses
    RewInc5P = 5;   %  five pulses (jackpot)

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
end

%% determine time interval before dimming should occur
% split the stimulus change/dimming times across conditions to allow for some control
num_cnd = length(unique(TrialRecord.ConditionsThisBlock)); % This variable is behaving strange and not as I understand its role. The use of unique ensures that it shows the actual number of conditions used in this block!
ccnd    = TrialRecord.CurrentCondition

itvtm   = linspace(dimmMin, dimmMax, num_cnd+1);
Dimmvec = dimmMin : dimmStep : dimmMax; % possible stimulus change times
if(num_cnd > 1)
    Dimmvec = Dimmvec(Dimmvec >= itvtm(ccnd) & Dimmvec < itvtm(ccnd+1)); % possible stimulus change times for current condition
end
cDimm   = Dimmvec(randi(length(Dimmvec)));  % stimulus change time used for the current trial

%  mu = (dimmMin + dimmMax)/2
%  max_range = exp(-(dimmMin/mu));
%  min_range = exp(-(dimmMax/mu));
%  R = dimmStep * round( (-mu*reallog((max_range-min_range)*rand(1,1)+min_range)) / dimmStep);

%% determine stimulus location
% location indices and labels
%  6  4  7   |   Tl  Tc  Tr
%  2  1  3   |   Ml  Mc  Mr
%  8  5  9   |   Bl  Bc  Br

if(~isfield(TrialRecord, 'LocVec'))
    TrialRecord.LocVec = [];   % track the joystick state across trials
end

if(FixGrid > 1 && ForceLoc == 0)
    if(isempty(TrialRecord.LocVec))
        TrialRecord.LocVec = repmat(1:FixGrid,1,LocRep);
    end

    TrialRecord.LocVec = TrialRecord.LocVec(randperm(length(TrialRecord.LocVec)));  % shuffle all possible locations

    cLoc = TrialRecord.LocVec(1);
elseif(ForceLoc > 0)
    cLoc = ForceLoc;
else
    cLoc = 1;
end

LocLbl = '--';
% get x coordinate
switch cLoc
    case {6,2,8}
        LocLbl(2) = 'l'
        StimX = -1 * FixStep;
    case {4,1,5}
        LocLbl(2) = 'c'
        StimX = 0;
    case {7,3,9}
        LocLbl(2) = 'r'
        StimX = FixStep;
end

% get y coordinate
switch cLoc
    case {6,4,7}
        LocLbl(1) = 'T'
        StimY = FixStep;
    case {2,1,3}
        LocLbl(1) = 'M'
        StimY = 0;
    case {8,5,9}
        LocLbl(1) = 'B'
        StimY = -1 * FixStep;
end

% move stimulus items to current location
reposition_object(StimGo,   StimX, StimY);
reposition_object(StimNoGo, StimX, StimY);

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

if(~isfield(TrialRecord, 'NoRespCnt'))
    TrialRecord.NoRespCnt = 0;  % track the trials with no behavioural response
    TrialRecord.ADDITI    = 0;  % Add extra wait time
end

disp('start timing file');

% check if lever remains released
[JoyHold] = eyejoytrack('holdtarget', JoyZero, JoyPullRad, 100); % this introduces an extra 100 ms interval

if(~JoyHold | TrialRecord.Jreleased == 0)
    disp('Wait for release');
    if(~JoyHold) % wait for release
        disp('is pressed');
        TrialRecord.Jreleased = 0;
        [JoyRel] = eyejoytrack('acquiretarget', JoyZero, JoyPullRad, pull_pause);

        if(JoyRel)
            eventmarker(39); % lever released
        end

    else % make sure that joystick is not pulled for a sufficient long time before starting a new trial
        disp('Wait for time out');
        [JoyHold] = eyejoytrack('holdtarget', JoyZero, JoyPullRad, pull_pause);

        if(JoyHold)
            TrialRecord.Jreleased = 1;
        else
            eventmarker(38); % lever pressed
        end
    end

    TrialRecord.CorrCount = 0;
    trialerror(8);
    cITI = 0;
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
    eventmarker(10);  % end ITI

    if(~isfield(TrialRecord, 'NextTrial'))
        TrialRecord.NextTrial = -1;
    end

%% #### Onset NoGo Stimulus #### %%
    disp('Start showing NoGo');

    NoGoOn = toggleobject(StimNoGo, 'EventMarker', [1, 35], 'Status', 'on');

    % get the computer time for the trial start as well
    NoGoEff = 1000 * toc(uint64(1));
    TrialRecord.Tstart(TrialRecord.CurrentTrialNumber) = NoGoEff;

    FixOn     = NoGoOn;
    TrialZero = NoGoOn;

%% #### Wait for Change #### %%
    disp(['Wait for dimming in ',int2str(cDimm),' ms']);
    [JoyHold PullTime] = eyejoytrack('holdtarget', JoyZero, JoyPullRad, cDimm);  % wait before stim change (hope no time passed since onset)
    EndWaitChng = 1000 * toc(uint64(1));

    % do not accept early responses
    if(~JoyHold)
        eventmarker(38);  % lever pressed

        % get the computer time for the lever release as well
        TPullEff = EndWaitChng;

        ErrCode  = 5;
        on_track = 0;
        StimOffTime = toggleobject(StimNoGo , 'EventMarker', 36,'Status', 'off');
        StimOffEff  = 1000 * toc(uint64(1));
        TrialRecord.CorrCount = 0;
        disp('Early release');
    end

%% #### Onset Go Stimulus #### %%
    % wait for release of joystick
    if(on_track)
        disp('Show Go Stimulus');
        GoOn = toggleobject([StimNoGo StimGo], 'EventMarker', [37]);  % show Go stimulus
        GoEff = 1000 * toc(uint64(1));

        [JoyHold PullTime] = eyejoytrack('holdtarget', JoyZero, JoyPullRad, max_resp);
        EndWaitResp = 1000 * toc(uint64(1));

        if(JoyHold)
            ErrCode  = 1; % no response
            on_track = 0;
            StimOffTime = toggleobject(StimGo, 'EventMarker', 36,'Status', 'off');
            StimOffEff  = 1000 * toc(uint64(1));
            TrialRecord.CorrCount = 0;

            % wait for possible late responses
            [JoyHold PullTime] = eyejoytrack('holdtarget', JoyZero, JoyPullRad, wait_late);
            EndWaitResp = 1000 * toc(uint64(1));

            if(~JoyHold)
                eventmarker(38);  % lever pressed
                ErrCode  = 2;     % late response
                TPullEff = EndWaitResp;
                RelTime  = TPullEff - GoEff;
                TrialRecord.NoRespCnt = 0;
            else
                TrialRecord.NoRespCnt = TrialRecord.NoRespCnt + 1;
            end
        else
            eventmarker(38);  % lever pressed

            TPullEff = EndWaitResp;
            RelTime  = TPullEff - GoEff;

            TrialRecord.NoRespCnt = 0;

            if(PullTime < wait_resp)
                ErrCode  = 5; % early response
                on_track = 0;
                StimOffTime = toggleobject(StimGo , 'EventMarker', 36,'Status', 'off');
                StimOffEff  = 1000 * toc(uint64(1));
            else
                rt = RelTime; % double check this one (use text table and compare to
            end
        end  % if(JoyHold)
    end  %  if(on_track)

%% #### Trial End #### %%
   eventmarker([2 17]);  % end of trial

    % finalize Trial
    if(on_track)  % correct trials
        ErrCode = 0;
        trialerror(0);  % correct

        TrialRecord.CorrCount = TrialRecord.CorrCount + 1;

        eventmarker(19);   % start pause
        idle(rew_lag);

        StimOffTime = toggleobject(StimGo, 'EventMarker', 36,'Status', 'off');
        StimOffEff = 1000 * toc(uint64(1));
        eventmarker([96 20]);

        RewEff = StimOffEff;
        goodmonkey(rew_dur, 'NumReward', rew_Npulse, 'PauseTime', rew_gap);
        eventmarker(96); % use this event twice to get an estimate of the reward duration, i.e. number of pulses.
    else   % error occurred
        trialerror(ErrCode);
        rew_Npulse = NaN;
        TrialRecord.CorrCount = 0;
        cITI = cITI + time_out;   % add a time-out to the current ITI
    end  %  if(on_track)

    FixOff    = StimOffTime - TrialZero;
    FixOffEff = StimOffEff  - NoGoEff;

    % delete current stimulus location from list if completed
    if(ErrCode ~= 8 && ErrCode ~= 9 && FixGrid > 1 && ForceLoc == 0)
        if(AccAll == 0 && on_track == 1)
            TrialRecord.LocVec(1) = [];
        elseif(AccAll == 1)
            TrialRecord.LocVec(1) = [];
        end
    end

    % get the computer time for the trial end as well
    TrialRecord.Tend(TrialRecord.CurrentTrialNumber) = 1000 * toc(uint64(1));

    % try to wait if lever is not yet released to avoid time out
    [JoyHold] = eyejoytrack('holdtarget', JoyZero, JoyPullRad, 100); % this introduces a 100 ms interval
    if(~JoyHold)
        disp('Wait release');

        [JoyRel] = eyejoytrack('acquiretarget', JoyZero, JoyPullRad, pull_pause);
        if(JoyRel)
            eventmarker(39); % lever released
            TrialRecord.Jreleased = 1;
        end
    end

%%%%###########################################%%%%
%%%%##########   Stop the Trial now   #########%%%%
%%%%###########################################%%%%
    % set ITI
    eventmarker([18 9]); % end post trial / start ITI

    if(TrialRecord.CurrentTrialNumber > 1)
        ITIeff = TrialRecord.Tstart(TrialRecord.CurrentTrialNumber) - TrialRecord.Tend(TrialRecord.CurrentTrialNumber-1);
    else
        ITIeff = NaN;
    end

    %% Determine ITI and check if the current Block needs to end with a pause
    if(~isfield(TrialRecord, 'BlockEnd'))
        TrialRecord.BlockEnd = toc(uint64(1)) + 60*BlockDur;
    end

    if(toc(uint64(1)) > TrialRecord.BlockEnd)
        TrialRecord.NoRespCnt = 0;
        breakvec = BlockBreakMin*60 : 10 : BlockBreakMax*60; % possible times for a break
        cbreak   = breakvec(randi(length(breakvec)));        % current break time

        cITI = cbreak*1000;
        TrialRecord.BlockEnd = toc(uint64(1)) + 60*BlockDur + cbreak;
        TrialRecord.NextTrial = TrialRecord.Tend(TrialRecord.CurrentTrialNumber) + cbreak*1000;
        disp(['Break Block for ',num2str(cbreak/60,'%.2f'), ' min']);
    else
        if(TrialRecord.NoRespCnt > 0 && ErrCode == 1)
            if(mod(TrialRecord.NoRespCnt,5) == 0 && TrialRecord.ADDITI < 45000 && useAddITI == 1)
                TrialRecord.ADDITI = TrialRecord.ADDITI + 5000;
            else
                TrialRecord.ADDITI = 0;
            end
            cITI = cITI + TrialRecord.ADDITI;
        else
            TrialRecord.NoRespCnt = 0;
            TrialRecord.ADDITI    = 0;
        end
        disp(['wait ITI ',int2str(cITI), ' ms']);
    end

    TrialRecord.NextTrial = TrialRecord.Tend(TrialRecord.CurrentTrialNumber) + cITI;
    set_iti(cITI);

    if(ErrCode ~= 8 && ErrCode ~= 9)
    %% write trial table
    % keep redundant information to check data and have some information as back up to reconstruct trials
    % Further, this text file should give fast and easy access to the behavioural performance. It is in a
    % format that allows to be read into R with the read.table command (use 'header=TRUE' option).

        tblnm = fullfile(tblpath, [cNHP,'_', cPRDGM,'_', cdt, '.dat']);

        if(~exist(tblpath,'dir'))
            mkdir(tblpath);
        end

        if(~exist(tblnm,'file'))  % create the table file
            tblptr = fopen(tblnm, 'w');

            fprintf(tblptr,'Date  Subject  Experiment  TrialNo  BlockNo  CondNo  StimLoc  StimX  StimY  StimSize  PreDimmGrey  PostDimmGrey  LocLbl  Result  TrialStart  TrialZero  FixOn  NoGoOn  GoOn  NoGoEff  GoEff  ChangeIntent  PullTime  TPullEff  RelTime  RT  StimOffTime  FixOff  FixOffEff  RewEff  NumRew  TrialEnd  ITIeff  NextITI  JoyPullRad\n');
        else
            tblptr = fopen(tblnm, 'a');
        end

        fprintf(tblptr, ...
            '%s  %s  %s  %6d  %4d  %4d  %d  %.4f  %.4f  %.4f  %.4f  %.4f  %s  %4d  %.4f  %6d  %6d  %6d  %6d  %.4f  %.4f  %6d  %6d  %.4f  %.4f  %.4f  %6d  %6d   %.4f  %.4f  %4d  %.4f  %.4f  %6d   %.4f\n', ...
            cdt, cNHP, cPRDGM, TrialRecord.CurrentTrialNumber, TrialRecord.CurrentBlock, ccnd, ...
            cLoc, StimX, StimY, sqrsz, NGsqrcol, Gsqrcol, LocLbl, ErrCode, ...
            TrialRecord.Tstart(TrialRecord.CurrentTrialNumber), TrialZero,FixOn, NoGoOn, GoOn, ...
            NoGoEff, GoEff, cDimm, PullTime, TPullEff, RelTime, rt, ...
            StimOffTime, FixOff, FixOffEff, RewEff, rew_Npulse, ...
            TrialRecord.Tend(TrialRecord.CurrentTrialNumber), ITIeff, cITI, JoyPullRad);

        fclose(tblptr);
    end
end  % if(on_track)  %% main

disp('end timing file');
