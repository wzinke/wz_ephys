function EV = TEMPO_EV_heitz
% interpretation done by wolf zinke.

%% define relevant channels
EV.EV_channel = 17;
EV.PD_channel =  2;

%% relevant event codes
EV.TrialStart_   = 1666;
EV.Target_       = 1666;  % this code seems to encode the trial start that apparently was defined by stimulus onset.
EV.Eot_          = 1667;
EV.TrialSuspend_ = 1668;

EV.Stimulation_  = 666;

EV.Abort         = 2621; % same as fixation break
EV.Abort_old     = 2620; % seems to be an old version to encode trial aborts

EV.FixSpotOn_  = 30; % only for D and E
EV.FixSpotOff_ = 40; % only for D and E
EV.Fixate_     = 31; % only for D and E

EV.Reward_     = 2727;
EV.Reward_tone = 2726;

%% task information
%This is encoded in duplets, first code defines type of
% informtion, subsequent code gives the value.
EV.Correct_     = 3034;
EV.Saccade_     = 3035; 
EV.TargetLoc    = 3013; 
EV.TargetCol    = 3008; 
EV.TargetClut   = 3009; % Min Target Color (CLUT value; higher values = BRIGHTER [more white])

EV.SetSize      = 3005; % code meaning: 0 -> set size = 2; 1 -> set size = 4; 2 -> set size = 8
EV.SetSizeCond  = 3004; % SET SIZE CONDITION (0 = 2; 1 = 4; 2 = 8; 3 = random (whatever this is supposed to mean)
EV.TaskType     = 3003; % 0 = fixation; 1 = detection; 2 = search/MG??
EV.HoldTime     = 3021; % use to determine search or memory guided (whatever this is supposed to mean)
EV.DistHomo     = 3032; % Homogeneous (0)/Non-homogeneous (1)
EV.Eccentricity = 3001; 

EV.FixTime      = 3026; % FixTime_Jit_ is the ms value monkey must be holding fixation (with some uniform or exponentially distributed jitter) before the target appears
EV.MG_Hold      = 3027; % Code actual MG hold time

% SAT specific task information
EV.SAT_Block      = 5000; % SAT blocking condtion (1 = blocked; 0 = random)
EV.SAT_Cond       = 5001; % SAT condition (1 = slow, 2 = med, 3 = fast)
EV.SAT_CutOff     = 5002; % SAT cutoff time for current block
EV.SAT_Percentile = 5003; % Percentile (may by unnecessary)
EV.SAT_Reward     = 5004; % Amount of reward for current block (time in ms solenoid on)
EV.SAT_Catch      = 3047; % Catch Trial %
EV.SAT_PunishLoc  = 6019; % Punish Time for Direction Error
EV.SAT_PunishTime = 6020; % Punish Time for missed deadline (noResp)
EV.SAT_SearchDur  = 3046; % Duration of SEARCH array (ms)
EV.SAT_FlipCond   = 6021; % Flip Conditions automatically every X trials

EV.SAT_DispClear  = 6022; % Was display cleared at deadline in FAST? if this code exist 1, otherwise 0

%% eye channel information
EV.GainsX    = 3036;
EV.GainsY    = 3037;
EV.Gains_XYF = 3038;
EV.Gains_YXF = 3039;

%% error codes
EV.FixError_        = 2621;
EV.CatchError       = 2622;
EV.HoldError        = 2623;
EV.LatencyError     = 2624;
EV.TargetHold_Error = 2625;
EV.Error_sacc       = 2626;
EV.Correct_too_fast = 2627;
EV.Correct_too_slow = 2628;

%% location specific codes
EV.SacLocVec  = 7000:7007;  % Expected codes for saccade locations
EV.StimLocVec = 6000:6007;  % code for information about the different stimuli (identifies location, but also if L or T and orientation)


%% get sorted arrays for event codes and event names
nmlst = fieldnames(EV);

for(i=1:length(nmlst))
    EV.names(1,i) = nmlst(i);
    EV.codes(1,i) = EV.(cell2mat(nmlst(i)))(1);
end

[EV.codes, idx] = sort(EV.codes);
EV.names = EV.names(idx);


