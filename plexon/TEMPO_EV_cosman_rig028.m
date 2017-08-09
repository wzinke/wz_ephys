function EV = TEMPO_EV_cosman_rig028

% define relevant channels
EV.EV_channel = 17;
EV.PD_channel =  2;

% relevant event codes

% Indicate the code block for trial infors
EV.StartInfos_ = 2998;
EV.EndInfos_   = 2999;
EV.InfosZero   = 3000;

% Trial start and trial end
EV.TrialStart_ = 1666;
EV.Eot_        = 1667;

% fixation spot and fixation state
EV.FixSpotOn_  = 2301;
EV.FixSpotOff_ = 2300;
EV.Fixate_     = 2660;

% Stimulus onset
EV.Target_     = 2651;

% response
EV.Saccade_    = 2810;
EV.Decide_     = 2811;

% feedback
EV.Reward_     = 2727;
EV.Reward_tone = 777;
EV.Tone_       = 2001;
EV.Error_tone  = 776;

% error codes
EV.Correct_    = 2600;
EV.FixError_   = 2750;
EV.GoError_    = 2752;
EV.Abort_      = 2620;

EV.GoTargFixError_    = 2754;
EV.PlacOn_            = 2320;
EV.StopSignal_        = 2653;
EV.GoSaccade_         = 2751;
EV.NOGOWrong_         = 2753;
EV.GOCorrect_         = 2755;
EV.NOGOCorrect_       = 2756;
EV.CatchCorrect_      = 2757;
EV.CatchIncorrectG_   = 2758;
EV.CatchIncorrectNG_  = 2759;
EV.BreakTFix_         = 2760;
EV.EarlySaccade_      = 2761;


EV.Error_sacc         = 887;
EV.Correct_sacc       = 888;

EV.Stimulation_       = 666;
EV.StimStart_         = 665;
EV.StimEnd_           = 667;


EV.ExtraReward_       = 2777;
EV.SoundOnReward_     = 2778;
EV.SoundNoReward_     = 2779;
EV.CmanHeader_        = 1501;
EV.MemHeader_         = 1502;
EV.GONOGOHeader_      = 1503;
EV.DelayedHeader_     = 1504;
EV.SearchHeader_      = 1507;
EV.Identify_Room_     = 1500;
EV.ZeroEyePosition_   = 2302;
EV.VSyncSynced_       = 999;
EV.MouthBegin_        = 2655;
EV.MouthEnd_          = 2656;
EV.MapHeader_         = 1503;
EV.FixWindow_         = 2770;
EV.TargetWindow_      = 2771;
EV.Staircase_         = 2772;
EV.Neg2Reinforcement_ = 2773;
EV.Feedback_          = 2774;
EV.RewardSize_        = 2927;
EV.TrialInBlock       = 2928;
EV.SendPenatrInfo_    = 2929;
EV.TargetPre_         = 2650;
EV.StopOn_            = 2654;
EV.StimFailed_        = 667;
EV.NPOS_              = 2721;
EV.Pos_               = 2722;
EV.SOUND_             = 2723;
EV.ISNOTNOGO_         = 2724;
EV.BIG_REWARD_        = 2725;
EV.TRIG_CHANGE_       = 2726;
EV.REWARD_RATIO_      = 2728;
EV.NOGO_RATIO_        = 2729;
EV.ISNOGO_            = 2730;
EV.STOP_ZAP_          = 2731;
EV.STIM_DUR_          = 2732;
EV.EXP_HOLDTIME_      = 2733;
EV.HOLDTIME_          = 2734;
EV.HOLD_JITTER_       = 2735;
EV.MAX_RESP_TIME_     = 2736;
EV.SOA_               = 2737;
EV.MAX_SOA_           = 2738;
EV.MIN_SOA_           = 2739;
EV.SOA_STEP_          = 2740;

%% get sorted arrays for event codes and event names
nmlst = fieldnames(EV);

for(i=1:length(nmlst))
    EV.names(1,i) = nmlst(i);
    EV.codes(1,i) = EV.(cell2mat(nmlst(i)));
end

[EV.codes, idx] = sort(EV.codes);
EV.names = EV.names(idx);


