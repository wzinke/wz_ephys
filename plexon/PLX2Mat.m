function plxMat = PLX2Mat(plxfile, EV, keep_sess, skip_waves, PreTrialTime, PostTrialTime, ITIstim, SampleRate)
% Read in a plexon data file and organize it according to trials.
%
% DESCRIPTION
%
%   This routine reads in a plexon data file and returns a Matlab struct that
%   contains the information in a trial sorted organization. Minimal knowledge
%   about the paradigm is used in order to keep this function more general. Task
%   information is kept with the event codes that could be used subsequently to
%   define the task outline and align trial times to a defined event.
%
%   To sort event codes and to define trial periods, an Event object is used that
%   has fields 'StartInfos_', 'EndInfos_', 'TrialStart_', and 'Eot_' that contain
%   the integer code that is associated with this event.
%
%   This routine needs the readPLXFileC function provided by B. Kraus on Matlab
%   File Exchange: http://www.mathworks.com/matlabcentral/fileexchange/42160-readplxfilec
%
%   The code is roughly based on the translation functions provided by R. Heitz and J. Cosman.
%
% SYNTAX
%
%   plxMat = PLX2Mat(plxfile, EV, keep_sess, skip_waves, PreTrialTime, PostTrialTime, ITIstim, SampleRate)
%
%   Input:
%
%       plxfile         - file name of the plexon data file
%
%       EV              - A structure that defines event codes or a tempo EVENTDEF file.
%                         This function works with a minimal set of event codes
%                         that must be present as field names in the EV struct:
%                             StartInfos_ - defines the begin of a trial info block
%                             EndInfos_   - defines the end of a trial info block
%                             TrialStart_ - defines the start of a trial
%                             Eot_        - defines the end of a trial
%
%       keep_sess       - keep analoge data and spike times as single vector for
%                         the complete session (memory demands increase considerably)
%                         > if set to 1 this outputs the session data and
%                           trial aligned data
%                         > if set to 2 it only outputs session data but
%                           does not do any trial sorting (no trial codes
%                           needed then!)
%
%       skip_waves      - do not keep single spike wave forms after the average
%                         is determined. This wil reduce memory usage a lot.
%
%       PreTrialTime    - keep data for this time period prior stimulus onset
%
%       PostTrialTime   - add this time to the end of a trial for prolonged sampling
%
%       ITIstim         - check for stimulation codes during the inter trial interval
%
%       SampleRate      - intended time reolution in [Hz]
%
% .........................................................................
% wolf zinke, wolfzinke@gmail.com
%
% $Created : 05-Jan-2015 by wolf zinke
%
% ToDo:  - Check the number of photo diode events and utilize the event they are
%          to achieve a more precise timing.
%
%        - consider using a filter/artefact detection on analog signal before trial splitting
%

% ____________________________________________________________________________ %
%% Set some default values
max_spont_rate = 150; % [spk/s] - maximum spontaneous activity to be accepted as valid DSP channel (avoid crappy noise channels)
min_resp_rate  =   6; % [spk/s] - minimum spontaneous and response activity to be accepted as valid DSP channel (avoid unresponsive units)
spontwin       = [-250   0];

% respwin        = [  50 200];
% min_trial = 40; % minimum number of trials to be processed

% ____________________________________________________________________________ %
%% Check if dependencies are fullfilled
if(~exist('readPLXFileC','file'))
    warning('readPLXFileC not found, please select now');
    [~,PathName] = uigetfile({'*.mex*'},'readPLXFileC mex file');
    addpath(PathName);
end

% ____________________________________________________________________________ %
%% get information about the current file version
% if(exist('wz_get_git_version','file'))
%     [gitversion, ~, GitInfo]    = wz_get_git_version(mfilename('fullpath'));
%     plxMat.TranslateFileVersion = gitversion;
%     plxMat.TranslateFileDate    = GitInfo.date;
% end

% ____________________________________________________________________________ %
%% check input data
% if not specified use GUI to get the file
if(~exist('plxfile','var') || isempty(plxfile))
    [FileName,PathName] = uigetfile({'*.plx;*.PLX'},'Load plexon file');
    plxfile = fullfile(PathName,FileName);
end

if(~exist('keep_sess','var') || isempty(keep_sess))
    keep_sess = 0; % [Hz] - keep times at this time resolution, but only use integer numbers for space reasons
end

if(~exist('skip_waves','var') || isempty(skip_waves))
    skip_waves = 0; % [Hz] - keep single spike waves
end

if(~exist('SampleRate','var') || isempty(SampleRate))
    SampleRate = 1000; % [Hz] - keep times at this time resolution, but only use integer numbers for space reasons
end

if(~exist('PreTrialTime','var') || isempty(PreTrialTime))
    PreTrialTime = -500; % [ms] - add this time before the actual trial start to sample data starting before the trial
end

if(~exist('PostTrialTime','var') || isempty(PostTrialTime))
    PostTrialTime = 100; % [ms] - add this time to the end of a trial for prolonged sampling
end

if(~exist('ITIstim','var') || isempty(ITIstim))
    ITIstim = 0; % flag to indicate whether stimulation occurred between trial (ad hoc hack for ultrasound experiments)
end

% ____________________________________________________________________________ %
%% define obj with key codes
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

% ____________________________________________________________________________ %
%% read in the plexon file
disp(['Reading plexon file: ',plxfile]); tic;
plx = readPLXFileC(plxfile,'all');
fprintf('... done after %.1f s\n\n', toc);

[plxMat.Path, plxMat.File] = fileparts(plxfile);

plxMat.DateNum = plx.Date;
plxMat.DateStr = datestr(plx.Date);

% ____________________________________________________________________________ %
%% determine event channels for event codes and photo diode strobes
if(~isfield(EV,'EV_channel'))
    cnt = 0;
    for(i=1:plx.NumEventChannels)
        if(length(plx.EventChannels(i).Timestamps) > 1)
            cnt = cnt +1;
            if(length(unique(plx.EventChannels(i).Values)) == 1)
                Pchan = i;
            else
                Tchan = i;
            end
        end
    end
    
    if(cnt > 2)
        warning('Something weired going on while determining event channels -> CHECK CAREFULLY!')
    end
else
    Tchan = EV.EV_channel;
    Pchan = EV.PD_channel;
end
% ____________________________________________________________________________ %
%% get information about the analog channels
% Identify active analog channels and get time of recording onset
plxMat.ActAD  = find([plx.ContinuousChannels(:).Enabled] == 1);
for(c=length(plxMat.ActAD):1) % make sure that active channels really contain data
    if(isempty(plx.ContinuousChannels(plxMat.ActAD(c)).Fragments))
        plxMat.ActAD(c) = [];
    end
end

plxMat.ADName     = {plx.ContinuousChannels(plxMat.ActAD).Name       };
plxMat.ADChannel  = [plx.ContinuousChannels(plxMat.ActAD).Channel    ];
plxMat.ADSmplRate = [plx.ContinuousChannels(plxMat.ActAD).ADFrequency];
plxMat.ADGain     = [plx.ContinuousChannels(plxMat.ActAD).ADGain     ];
plxMat.ADPreGain  = [plx.ContinuousChannels(plxMat.ActAD).PreAmpGain ];
plxMat.ADStartRec = 1000 * [plx.ContinuousChannels(plxMat.ActAD).Timestamps] / plx.ADFrequency; % nout converting to double is a trick to round to full ms

% This is IMPORTANT information: The start of sampling analog data does not
% define time zero in plexon. Instead, this start occurs with a lag to the
% internal zero time that is used to time stamp event codes and spike times. The
% analog channels also have a time stamp that indicates the start of sampling
% relative to the event zero time. If not accounted for this the analog signal
% will be shifted relative to event times and spike times.
% ADRecOn = unique([plx.ContinuousChannels(plxMat.ActAD).Timestamps]) / plx.ADFrequency;

if(length(unique(plxMat.ADSmplRate)) > 1)
    error('Sampling rate of analog channels differs! This needs an adjustment of the code!');
end

plxMat.SessTimeMinutes = double(plx.LastTimestamp) / (60*plx.ADFrequency);

% ____________________________________________________________________________ %
%% get event timings
% convert to timing in milliseconds and retain double precision;
% subtract ADRecOn from the event timings to ensure that times are aligned with
% the first time bin of the analog signal.
plxMat.PD_Times    = SampleRate * double(plx.EventChannels(Pchan).Timestamps) / plx.ADFrequency; % photo diode
plxMat.Event_Times = SampleRate * double(plx.EventChannels(Tchan).Timestamps) / plx.ADFrequency; % event time stamps
plxMat.Event_Codes = int16(plx.EventChannels(Tchan).Values);  % codes identifying the events corresponding to the time stamps

% initcodes  = plxMat.Event_Codes(1:10);
% plxMat.rig = initcodes(initcodes<100);% hope ths was used consistently to encode room number

% Get available event codes
EVlst = sort(unique(plxMat.Event_Codes));
[~,ia,ib] = intersect(EVlst, EV.codes);

plxMat.EV.codes = int16(EV.codes(ib));
plxMat.EV.names = EV.names(ib);

plxMat.EV.unspec = int16(EVlst);
plxMat.EV.unspec(ia) = [];

for(i=1:length(plxMat.EV.codes))
    plxMat.EV.(plxMat.EV.names{i}) = plxMat.EV.codes(i);
end

% ____________________________________________________________________________ %
%% re-order Event codes and times into a trial matrix
if(keep_sess ~= 2)
    disp(['Determining event trial positions ']); tic;
    
    % get trial event codes
    if(isfield(plxMat.EV,'TrialStart_'))
        Tstart = find(plxMat.Event_Codes == plxMat.EV.TrialStart_);  % identifies all trial chunk initialized by tempo
        Tend   = find(plxMat.Event_Codes == plxMat.EV.Eot_);
    else
        Tstart = [];
    end
    
    if(isempty(Tstart))
        RestData  = 1;  % resting state/movie data with no trials - treat whole session as single tria
        keep_sess = 0; % no reason to keep session if it is the single trial anyway.
        Tstart    = 1;
        Tend = length(plxMat.Event_Codes);
    else
        RestData = 0;
        
        % check that trial events are complete and not chopped at start or end
        if(length(Tstart) == length(Tend)+1 && Tstart(end) > Tend(end))
            warning('No trial end code found for last trial. Will discard last trial!');
            Tstart(end) = [];
        end
        
        if(length(Tstart)+1 == length(Tend) && Tstart(1) > Tend(1))
            warning('No trial start code found for first trial. Will discard first trial end!');
            Tend(1) = [];
        end
        
        % if there is a mismatch in number of codes correct this error
        if(length(Tstart) > length(Tend))
            warning('Apparently, trial end codes were ommited during the session. Try to identify the trials and remove it.');
            tmpTS  = nan(length(Tend),1);
            for(ct = 1:length(Tend))
                p = find(Tstart < Tend(ct),1,'last');
                tmpTS(ct) = Tstart(p);
            end
            Tstart = tmpTS;
            clear tmpTS;
        end
        
        if(length(Tstart) < length(Tend))
            warning('Apparently, trial start codes were ommited during the session. Try to identify the trials and remove it.');
            tmpTE = nan(length(Tstart),1);
            for(ct = 1:length(Tstart))
                p = find(Tend > Tstart(ct),1,'first');
                tmpTE(ct) = Tend(p);
            end
            Tend = tmpTE;
            clear tmpTE;
        end
        
        p = find(diff(Tend) == 0);
        if(~isempty(p))
            Tend(p-1)   = [];
            Tstart(p-1) = [];
        end
        
        p = diff(Tstart) == 0;
        if(sum(p) > 0)
            Tend(p)   = [];
            Tstart(p) = [];
        end
    end
else
    RestData = 1;  % resting state/movie data with no trials - treat whole session as single tria
    Tstart   = 1;
    Tend  = length(plxMat.Event_Codes);
end

% % this loop is required because some data files have for some reason more end of
% % trials (Eot_) encoded than trials starts (TrialStart_).
% Tend = nan(length(Tstart),1);
% for(tn = 1:length(Tstart))
%     tep = tmpTend(find(tmpTend > Tstart(tn),1, 'first'));
% %     tep = find(plxMat.Event_Codes == plxMat.EV.Eot_ & plxMat.Event_Times > plxMat.Event_Times(Tstart(tn)),1,'first');
%     if(~isempty(tep))
%         Tend(tn) = tep;
%     end
% end

if(RestData == 0)
    % get trial information codes
    if(isfield(plxMat.EV,'StartInfos_')) % blocks with trial related information are clearly specified
        % identifies all info chunks initialized by tempo
        tmpIstart = find(plxMat.Event_Codes == plxMat.EV.StartInfos_) + 1;
        tmpIend   = find(plxMat.Event_Codes == plxMat.EV.EndInfos_)   - 1;
        
        % pre-allocate
        Istart = nan(length(Tstart), 1);
        Iend   = nan(length(Tstart), 1);
        rem_trial=[];
        
        % loop over all trial starts to identify subsequent info blocks for each trial
        for(ct = 1:length(Tstart))
            if(ct < length(Tstart))
                Ps = find(tmpIstart > Tend(ct) & tmpIstart < Tstart(ct+1), 1, 'first');
                Pe = find(tmpIend   > Tend(ct) & tmpIend   < Tstart(ct+1), 1, 'first');
            else
                Ps = find(tmpIstart > Tend(ct), 1, 'first');
                Pe = find(tmpIend   > Ps,       1, 'first');
            end
            
            if(isempty(Ps) || isempty(Pe))
                rem_trial = [rem_trial; ct];
            else
                Istart(ct) = tmpIstart(Ps);
                Iend(ct)   = tmpIend(Pe);
            end
        end
        
        clear tmpIstart tmpIend
        
        % delete trials if no info block was found
        if(~isempty(rem_trial))
            warning('Not all trials had an info block! Removed following trials:');
            disp(rem_trial);
            Tstart(rem_trial) = [];
            Tend(rem_trial)   = [];
            Istart(rem_trial) = [];
            Iend(rem_trial)   = [];
        end
    else % no clean separation into event codes and info codes (old T/L search data and SAT data).
        % Some trial relevant informations are stored prior the trial, during the trial,
        % but mainly afterwards. To keep most of the codes, the info block
        % contains everything from the trial start to some arbitrary trial
        % end or info end event. This will be dirty...
        Istart = Tstart;
        Iend   = nan(size(Istart));
        
        % Determine end of trial chunk. This is really ugly...
        for(tn = 1:length(Iend))
            p1 = 1:length(plxMat.Event_Times) > Tend(tn)+1;
            p2 = plxMat.Event_Codes == plxMat.EV.TrialStart_ | plxMat.Event_Codes == plxMat.EV.Eot_; % check later if there are any trial abort codes in this chunk!
            
            pe = find(p1' == 1 & p2 == 1 , 1, 'first') - 1 ;
            
            if(~isempty(pe))
                Iend(tn) = pe;
            else
                if(tn == length(Istart))
                    Iend(tn) = length(plxMat.Event_Codes);
                end
            end
            % BE AWARE: this trial chunk is not accurate: there are codes
            % preceding the trial start code that belong to the subsequent trial.
            % These old data files do not contain clear codes that identify
            % accurately a complete trial block!
        end
    end
end

if(RestData == 0)
    % get photo diode refresh times for each trial
    PDpos = [bsxfun(@ge, plxMat.PD_Times(:)', plxMat.Event_Times(Tstart(:))) & ...
        bsxfun(@le, plxMat.PD_Times(:)', plxMat.Event_Times(Tend(:)))];
    
    % a valid trial must contain a stimulus onset as defined in the second column  of the
    % PDpos variable. Identify trials without stimulus presentation and remove them.
    pdcnt = sum(PDpos,2);
    %     rem_invalid = find(pdcnt<2);
    %
    %     if(~isempty(rem_invalid))
    %
    %         Istart(rem_invalid,:) = [];
    %         Iend(  rem_invalid,:) = [];
    %
    %         Tstart(rem_invalid,:) = [];
    %         Tend(  rem_invalid,:) = [];
    %         PDpos( rem_invalid,:) = [];
    %         pdcnt = sum(PDpos,2);
    %     end
    %
    % derive more information from the event blocks
    NITrials = length(Istart); % define a trial as chunk with stimulus presentation
    
    %     if(NITrials < min_trial)
    %         plxMat = [];
    %         warning('Not sufficient trials! --> skipped.');
    %     end
    
    % get event codes for trial information associations with trials
    infocnt = Iend - Istart + 1;
end

evcnt = Tend - Tstart + 1;
NTrialEV = max(evcnt);

plxMat.NTrials = length(Tstart);

fprintf('... done after %.1f s\n\n', toc);

% ____________________________________________________________________________ %
%% sanity checks
if(RestData == 0)
    if(NITrials ~= plxMat.NTrials)
        if(abs(NITrials - plxMat.NTrials) == 1)
            if(Tstart(1) > Istart(1))
                Istart(1) = [];
                Iend(1)   = [];
                warning('First trial info recorded before first trial start code! - first trial skipped.');
            elseif(Tstart(end) > Istart(end))
                Tstart(end) = [];
                Tend(end)   = [];
                warning('Last trial without trial info! - last trial skipped.');
            else
                error('Unexpectected mismatch of number of trials. -> CHECK CAREFULLY!')
            end
        else
            error('Something weird going on: different number of info blocks and trials -> CHECK CAREFULLY!')
        end
    end
    
    
    if(any(plxMat.Event_Times(Tstart) > plxMat.Event_Times(Istart)))
        error('Something weird going on: Info block found before trial start -> CHECK CAREFULLY!')
    end
    
    if(isempty(Tstart))
        warning('No trial start found -> CHECK CAREFULLY!')
    elseif(length(Tstart) ~= length(Tend))
        error('Unequal number of info trial start codes and trial end codes -> CHECK CAREFULLY!')
    end
    
    if(sum(pdcnt) == 0)
        warning('No photo diode encodes found -> CHECK CAREFULLY!')
    elseif(size(PDpos,1) ~= plxMat.NTrials)
        warning('Not all trials have a photo diode event -> CHECK CAREFULLY!')
    end
end

% ____________________________________________________________________________ %
%% create a trial matrix (rows correspond to trials)
disp('creating Event trial matrices '); tic;

if(RestData == 0)
    plxMat.TrialMat_INFOcodes = nan(NITrials, max(infocnt));
    
    plxMat.TrialMat_EVcodes = nan(plxMat.NTrials, NTrialEV);
    plxMat.TrialMat_EVtimes = nan(plxMat.NTrials, NTrialEV);
    plxMat.TrialMat_PDtimes = nan(plxMat.NTrials, max(pdcnt));
    
    for(t = 1:plxMat.NTrials)
        plxMat.TrialMat_INFOcodes(t, 1:infocnt(t)) = plxMat.Event_Codes(Istart(t):Iend(t))';
        plxMat.TrialMat_EVcodes(  t, 1:evcnt(t))   = plxMat.Event_Codes(Tstart(t):Tend(t))';
        plxMat.TrialMat_EVtimes(  t, 1:evcnt(t))   = plxMat.Event_Times(Tstart(t):Tend(t))';
        if(pdcnt(t) > 0)
            plxMat.TrialMat_PDtimes(t, 1:pdcnt(t)) = plxMat.PD_Times(PDpos(t,:));
        end
    end
end

if(keep_sess ~= 2)
    % get the time of the event that defines time zero of the trial
    % TrialZero = plxMat.TrialMat_PDtimes(:,2); % this was supposed to be stimulus onset
    TrialZero = plxMat.Event_Times(Tstart); % allign to trial onset
    
    % get trial times relative to first trial start (diff(plxMat.SessTime) will give inter-stimulus interval)
    plxMat.SessTime = TrialZero - TrialZero(1);
    
    % get trial start and trial end as session time
    plxMat.TrialStartEndDur      = [TrialZero, plxMat.Event_Times(Tend)];
    plxMat.TrialStartEndDur(:,3) = diff(plxMat.TrialStartEndDur,1,2);
    
    % allign trials to the event that defines time zero
    plxMat.TrialMat_PDtimes = bsxfun(@minus, plxMat.TrialMat_PDtimes, TrialZero);
    plxMat.TrialMat_EVtimes = bsxfun(@minus, plxMat.TrialMat_EVtimes, TrialZero);
    
    fprintf('... done after %.1f s\n\n', toc);
    
    if(keep_sess == 0) % clean up some data that is not wanted
        if(ITIstim == 0)
            plxMat = rmfield(plxMat, 'Event_Codes');
            plxMat = rmfield(plxMat, 'Event_Times');
        end
        plxMat = rmfield(plxMat, 'PD_Times');
    end
end
clear infopos trialpos PDpos infocnt evcnt pdcnt % just reduce the memory load a bit

% ____________________________________________________________________________ %
%% Analog data
disp(['Getting trial matrices of analog data']); tic;

if(keep_sess ~= 2)
    Bstart = round(plxMat.TrialStartEndDur(:,1)) + PreTrialTime;
    Bend   = round(plxMat.TrialStartEndDur(:,2)) + PostTrialTime;
    
    bincnt         = Bend-Bstart+1;
    NTrialBin      = max(bincnt);
    plxMat.timevec = [0:NTrialBin-1] + PreTrialTime; % aligned to TrialZero
end

% this will eat time: loop over channels and trials
for(c=1:length(plxMat.ActAD))
    cchan = sprintf('AD%02d', plxMat.ActAD(c));
%     cchan = cell2mat(plxMat.ADName(c));
    
    % create a full time vector. This should make sure that each entry in the verctor
    % for the analog data corresponds to a 1 ms step. Initial values are filled with
    % NaNs to compensate for the lag when recording started. If pauses occured during
    % recording, the time between two blocks also is filled with NaN values.
    
    TimeStamps = round(SampleRate * double(plx.ContinuousChannels(plxMat.ActAD(c)).Timestamps)  / plx.ADFrequency);
    
    cADval = nan(1,TimeStamps(end) + plx.ContinuousChannels(plxMat.ActAD(c)).Fragments(end));
    
    if(length(cADval) < Bend(end)) % in case there was not sufficient time between last trial and end of daq.
        cADval(length(cADval)+1 : Bend(end)) = NaN;
    end
    
    ADstart = 1;
    for(a = 1:length(TimeStamps))
        cST  = TimeStamps(a);
        if(cST == 0)
            cST = 1;
        end
        cADN = plx.ContinuousChannels(plxMat.ActAD(c)).Fragments(a);
        cADval(cST:cST+cADN-1) = double(plx.ContinuousChannels(plxMat.ActAD(c)).Values(ADstart:ADstart+cADN-1));
        ADstart = sum(plx.ContinuousChannels(plxMat.ActAD(c)).Fragments(1:a)) + 1;
    end
    
    if(keep_sess ~= 0)
        plxMat.([cchan,'_session']) = cADval;
    end
    
    %     % apply z-scoring
    %     AD_z = zscore(double(plx.ContinuousChannels(plxMat.ActAD(c)).Values));
    
    if(keep_sess ~= 2)
        plxMat.(cchan) = nan(plxMat.NTrials, NTrialBin);
        
        for(t=1:plxMat.NTrials)
            plxMat.(cchan)(t,1:bincnt(t)) = cADval(Bstart(t):Bend(t));
            %         plxMat.([cchan,'_z'])(t,1:bincnt(t)) = AD_z(Bstart(t):Bend(t));
        end
        
        if(c==1)
            plxMat.trialcnt_vec = sum(~isnan(plxMat.(cchan)));
        end
    end
end
fprintf('... done after %.1f s\n\n', toc);

clear timebins bincnt cADval % just reduce the memory load a bit

% ____________________________________________________________________________ %
%% Get spike times
disp(['Getting trial matrices of spike times']); tic;
spk_sfx = 'a':'k'; % there should not be much more needed

% identify channels that contain unit activity
unitcnt = sum(plx.SpikeTimestampCounts > 0);
spkchan = find(unitcnt > 0);

plxMat.SpikeChan = spkchan;

% assignment of electrode type; This is a bit ad hoc and should primarily help
% for sanity checks. Given the number of electrode channels it might be guessed
% what electrode was used. However, this also depends on the assignment of the
% channels to the elctrode.
if(max(spkchan) > 24)
    plxMat.ElectrodeType = 'neuronexus';
elseif(max(spkchan) > 10)
    plxMat.ElectrodeType = 'Uprobe';
elseif(length(spkchan)  == 1 )
    plxMat.ElectrodeType = 'Single';
elseif(max(spkchan) < 8)
    plxMat.ElectrodeType = 'Array';
else
    warning('Could not figure out what electrode type was used!');
    plxMat.ElectrodeType = 'NA';
end

plxMat.DSP_names = {};
plxMat.num_units =  0;

for(s=1:length(spkchan))
    unitlst = unique(plx.SpikeChannels(spkchan(s)).Units);
    unitlst(unitlst==0) = [];  % 0 = unsorted unit, 1,2,3,4, ... = sorted units a,b,c,d, ... - just get rid of unsorted spikes
    
    if(~isempty(unitlst))
        cspiketimes = SampleRate * double(plx.SpikeChannels(spkchan(s)).Timestamps) / plx.ADFrequency;
        
        if(skip_waves == 0)
            unitwaves = double(plx.SpikeChannels(spkchan(s)).Waves)'; % get waveform matrix
            norm_fac  = median(abs(unitwaves(:))); % assuming that the mean signal should be zero
        end
        
        for(u=1:length(unitlst))
            spkpos = plx.SpikeChannels(spkchan(s)).Units == unitlst(u);
            unitspikes = cspiketimes(spkpos);
            
            cunit = sprintf('DSP%02d%c', spkchan(s), spk_sfx(unitlst(u)));
            if(keep_sess ~= 0)
                plxMat.DSP.([cunit, '_session']) = unitspikes;
            end
            
            if(keep_sess ~= 2)
                spkcnt = zeros(plxMat.NTrials,1);
                spktms = [];
                trial_spikes = [];  % select only spike waveforms occurring within a valid trial - reduce file size and keep only relevant data.
                spike_trials = [];
                
                for(t=1:plxMat.NTrials)
                    pos = find(unitspikes >= plxMat.TrialStartEndDur(t,1) & ...
                        unitspikes <= plxMat.TrialStartEndDur(t,2)+PostTrialTime);
                    
                    if(~isempty(pos))
                        spkcnt(t,1) = length(pos);
                        spktms(t, 1:spkcnt(t)) = unitspikes(pos);
                        trial_spikes = [trial_spikes; pos(:)];
                        spike_trials = [spike_trials ; repmat(t,spkcnt(t),1)];
                    else
                        spktms(t, 1:spkcnt(t)) = NaN;
                    end
                end
                
                % matlab fills the matrix with zeroes, get rid of them here
                spktms(spktms == 0) = NaN;
                
                % align to stimulus onset. Do it after removing 0 values to not
                % remove spikes that indeed  did occurr exactly on stimulus onset
                spktms = bsxfun(@minus, spktms, TrialZero);
                
                trialrate = SampleRate .* spkcnt./plxMat.TrialStartEndDur(:,3); % get the average firing rate per trial
                mediantrialrate = median(trialrate(trialrate > 0));
                if(isempty(mediantrialrate))
                    mediantrialrate = 0;
                end
                
                %             respspk   = sum(spktms(:) >= respwin(1) & spktms(:) <= respwin(2));
                %             respspkrt = 1000 * respspk / plxMat.NTrials / (diff(respwin)+1);
                
                spospk   = sum(spktms(:) >= spontwin(1) & spktms(:) <= spontwin(2));
                spospkrt = 1000 * spospk / plxMat.NTrials / (diff(spontwin)+1);
                
                if(spospkrt < max_spont_rate && mediantrialrate > min_resp_rate && max(spkcnt) > 4)   %&& mediantrialrate > min_max_rate(1) && mediantrialrate < min_max_rate(2)) % discard empty channels and crap channels
                    plxMat.num_units = plxMat.num_units + 1;
                    plxMat.DSP_names(plxMat.num_units,1) = {cunit};
                    
                    plxMat.DSP.(cunit) = spktms;
                    
                    if(skip_waves == 0)
                        plxMat.DSP.([cunit, '_wave']).unit         = cunit;
                        plxMat.DSP.([cunit, '_wave']).channel      = plx.SpikeChannels(spkchan(s)).Channel;
                        plxMat.DSP.([cunit, '_wave']).suffix       = spk_sfx(unitlst(u));
                        plxMat.DSP.([cunit, '_wave']).thresh       = double(plx.SpikeChannels(spkchan(s)).Threshold) ./ norm_fac;
                        plxMat.DSP.([cunit, '_wave']).norm_fac     = norm_fac;
                        plxMat.DSP.([cunit, '_wave']).num          = length(trial_spikes);
                        plxMat.DSP.([cunit, '_wave']).spike_trials = spike_trials;
                        plxMat.DSP.([cunit, '_wave']).waves        = unitwaves(trial_spikes,:) ./ norm_fac;
                    end
                else
                    if(keep_sess == 1)
                        plxMat.DSP = rmfield(plxMat.DSP,([cunit, '_session']));
                    end
                end
            end
        end  % for(u=1:length(unitlst))
    end  % if(~isempty(unitlst))
end  % for(s=1:length(spkchan))

fprintf('... done after %.1f s\n\n', toc);

% ____________________________________________________________________________ %
%%
if(exist('ByetSize','file'))
    fprintf('Resulting structure has a size of %s', ByetSize(plxMat));
end
