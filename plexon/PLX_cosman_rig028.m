% script PLX_cosman_rig028

%% predefine expected Task variables
plxMat.Task.TargetLoc     = nan(plxMat.NTrials,1);  % Target location (angle)
plxMat.Task.TargetPos     = nan(plxMat.NTrials,1);  % Target location (index number)
plxMat.Task.SOA           = nan(plxMat.NTrials,1);  % ?

if(any(plxMat.TrialMat_INFOcodes(:) == 3000))  % MG trials exist
    plxMat.Task.MGcolor       = nan(plxMat.NTrials,1);  % color of item in the MG task
    plxMat.Task.DelayDur      = nan(plxMat.NTrials,1);  % Delay duration for memory task
end

if(any(plxMat.TrialMat_INFOcodes(:) == 3999))  % Search task (including capture) trials exist
    plxMat.Task.ArrayStruct   = nan(plxMat.NTrials,1);  % Structured vs unstrcutured display (typical vs. contextual cue)
    plxMat.Task.DistHomo      = nan(plxMat.NTrials,1);  % homo/hetero/random search cond
    plxMat.Task.SingMode      = nan(plxMat.NTrials,1);  % Could a singleton be expected in a trial?
    plxMat.Task.SetSize       = nan(plxMat.NTrials,1);  % Number of items in search array
    plxMat.Task.TargetType    = nan(plxMat.NTrials,1);  % Target Type - T (1) vs L (2)
    plxMat.Task.TrialType     = nan(plxMat.NTrials,1);  % Trial Type (random vs repeated displays)
    plxMat.Task.Eccentricity  = nan(plxMat.NTrials,1);  % Eccentricity of search array items
    plxMat.Task.Singleton     = nan(plxMat.NTrials,1);  % singleton presence: 1111 - singleton absent; 2222 - singleton present
    plxMat.Task.TargetHemi    = nan(plxMat.NTrials,1);  % Target hemifield: 8200 - right; 8100 - left;8888 - midline
    plxMat.Task.DistHemi      = nan(plxMat.NTrials,1);  % Distractor hemifield: 8200 - right; 8100 - left;8888 - midline
    plxMat.Task.DistLoc       = nan(plxMat.NTrials,1);  % Distractor location (angle)
    plxMat.Task.DistPos       = nan(plxMat.NTrials,1);  % Distractor location (index number)
    plxMat.Task.IsCatch       = nan(plxMat.NTrials,1);  % current trial is catch trial
    plxMat.Task.SingleDistCol = nan(plxMat.NTrials,1);  % Singleton Distractor Color
    plxMat.Task.DistOri       = nan(plxMat.NTrials,1);  % Distractor Orientation: 0-3 (0 = upright, 1= right tilt, 3 = inverted, 4 = left tilt)
    plxMat.Task.TargetOri     = nan(plxMat.NTrials,1);  % Target Orientation:     0-3 (0 = upright, 1= right tilt, 3 = inverted, 4 = left tilt)
    plxMat.Task.TargetHold    = nan(plxMat.NTrials,1);  % time to keep gaze at target to be accepted as correct response
    plxMat.Task.PercSgl       = nan(plxMat.NTrials,1);  % Percent trials w/ singleton distractor
    plxMat.Task.PercCatch     = nan(plxMat.NTrials,1);  % Percent catch trials
end

%% loop over trials and identify relevant information
for(t=1:plxMat.NTrials)
    cInf = plxMat.TrialMat_INFOcodes(t,:);  % get the info codes corresponding to the current trial

    if(cInf(1) == 3000)
        plxMat.Task.TaskType{t} = 'MG';  % memory guided saccade task
    elseif(cInf(1) == 3999)
        plxMat.Task.TaskType{t} = 'Search';
    else
        plxMat.Task.TaskType{t} = 'NA';
    end

    % the numerical values that are subtracted are defined in the tempo
    % configuration files. Unfortunately this has to be hard coded and is
    % arbitrary in its definition. Be careful and thouroughy double check.
    switch char(plxMat.Task.TaskType(t))
        case 'MG'
            plxMat.Task.TargetLoc(t) = cInf(19) - 3000;
            plxMat.Task.MGcolor(t)   = cInf(20) - 3000;

            % get the delay duration based on fixation spot offset
            plxMat.Task.DelayDur(t) = plxMat.Task.FixSpotOff(t) - plxMat.Task.StimOnset(t);

        case 'Search'
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

            if(cInf(10) == 8888)
                plxMat.Task.TargetHemi(t) = 'm';
            elseif(cInf(10) == 8200)
                plxMat.Task.TargetHemi(t) = 'r';
            elseif(cInf(10) == 8100)
                plxMat.Task.TargetHemi(t) = 'l';
            end

            if(cInf(11) == 8888)
                plxMat.Task.DistHemi(t) = 'm';
            elseif(cInf(11) == 8200)
                plxMat.Task.DistHemi(t) = 'r';
            elseif(cInf(11) == 8100)
                plxMat.Task.DistHemi(t) = 'l';
            end

            if(plxMat.Task.SetSize(t) == 1) % with only one item it is a detection task
                plxMat.Task.TaskType{t} = 'Det';  % detection task
            elseif(plxMat.Task.SingMode(t) == 1)  % double check this one
                plxMat.Task.TaskType{t} = 'Cap';  % capture task with salient singleton distractor
            end

        plxMat.Task.TargetPos(t) = mod((360-plxMat.Task.TargetLoc(t)+90)/45,8); % transform angle to index position
        plxMat.Task.DistPos(t)   = mod((360-plxMat.Task.DistLoc(t)  +90)/45,8);
    end

    plxMat.Task.MaxSaccDur(t)   = cInf(23) - 3000;
    plxMat.Task.MaxSaccTime(t)  = cInf(24) - 3000;
    plxMat.Task.TimeOut(t)      = cInf(25) - 3000;
    plxMat.Task.RewardDur(t)    = cInf(26) - 3000;
    plxMat.Task.RewardOffset(t) = cInf(27) - 3000;
    plxMat.Task.TargetHold(t)   = cInf(28) - 3000;
    plxMat.Task.Gains_EyeX(t)   = cInf(30) - 1000;
    plxMat.Task.Gains_EyeY(t)   = cInf(31) - 1000;
    plxMat.Task.Gains_XYF(t)    = cInf(32) - 1000;
    plxMat.Task.Gains_YXF(t)    = cInf(33) - 1000;
end
