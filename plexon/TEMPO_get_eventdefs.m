function [EVCODE] = TEMPO_get_eventdefs(EVDEFFL)
% Get Event codes from a tempo event definition file
%  
% DESCRIPTION 
%
% SYNTAX 
% 
%   [EVCODE] = TEMPO_get_eventdefs(EVDEFFL)
%
%   Input:
%
%       EVDEFFL  File name of the event definition file. If not specified, a GUI
%                will be opened that allows to manually select this file.
%
%   Output:
%
%       EVCODE   A matlab struct that contains two arrays, one vector of integer
%                numbers and a corresponding cell array that specify the meaning
%                of these numbers. In addition, all codes are stored as seperate
%                fields in the struct to allow using meaningful names instead of
%                crypric numbers.
%
% ......................................................................... 
% wolf zinke, wolfzinke@gmail.com 
%
% $Created : 22-Jan-2015 by wolf zinke

if(~exist('EVDEFFL','var') || isempty(EVDEFFL))
    [FileName,PathName] = uigetfile({'*.pro'},'tempo event definition file');
    EVDEFFL = fullfile(PathName,FileName);
end

if(~exist(EVDEFFL,'file'))
    error('File not found!');
end

% i am too lacy to implement this in Matlab ...
system(['cat ',EVDEFFL, ' | grep "declare"  | cut -d";" -f1 | sed -e "s/declare hide constant //g" |', ...
        'sed -e "s/declare hide float //g" | grep -v // | sed -e "s/[[:space:]]\+//g" | ' ...  '
        'awk "{\$1=\$1}{ print }" | sed -e "s/ //g" | sed -e "s/=/ /g" > tmp_ev.dat']);   
    
fid = fopen('tmp_ev.dat');
defs = textscan(fid,'%s %d');
fclose(fid); delete('tmp_ev.dat');

EVCODE.File  = EVDEFFL;
EVCODE.names = defs{1};
EVCODE.codes = double(defs{2});

[EVCODE.codes, idx] = sort(EVCODE.codes);
EVCODE.names = EVCODE.names(idx);

for(i=1:length(EVCODE.names))
    EVCODE.(EVCODE.names{i}) = EVCODE.codes(i);
end

%% define error codes
% TODO TODO!!!
% define a extra field for error codes and names. Look for possible code
% entries, name them in an interpretable way and store name and code in order to
% identify later these existing errors.
% TODO TODO!!!

%% Example code set
%            Stimulation_: 666
%             StimFailed_: 667
%              Error_tone: 776
%             Reward_tone: 777
%              Error_sacc: 887
%            Correct_sacc: 888
%            VSyncSynced_: 999
%          Identify_Room_: 1500
%             CmanHeader_: 1501
%              MemHeader_: 1502
%           GONOGOHeader_: 1503
%              MapHeader_: 1503
%          DelayedHeader_: 1504
%           SearchHeader_: 1507
%             TrialStart_: 1666
%                    Eot_: 1667
%                   Tone_: 2001
%             FixSpotOff_: 2300
%              FixSpotOn_: 2301
%        ZeroEyePosition_: 2302
%                 PlacOn_: 2320
%                Correct_: 2600
%                  Abort_: 2620
%              TargetPre_: 2650
%                 Target_: 2651
%             StopSignal_: 2653
%                 StopOn_: 2654
%             MouthBegin_: 2655
%               MouthEnd_: 2656
%                 Fixate_: 2660
%                   NPOS_: 2721
%                    Pos_: 2722
%                  SOUND_: 2723
%              ISNOTNOGO_: 2724
%             BIG_REWARD_: 2725
%            TRIG_CHANGE_: 2726
%                 Reward_: 2727
%           REWARD_RATIO_: 2728
%             NOGO_RATIO_: 2729
%                 ISNOGO_: 2730
%               STOP_ZAP_: 2731
%               STIM_DUR_: 2732
%           EXP_HOLDTIME_: 2733
%               HOLDTIME_: 2734
%            HOLD_JITTER_: 2735
%          MAX_RESP_TIME_: 2736
%                    SOA_: 2737
%                MAX_SOA_: 2738
%                MIN_SOA_: 2739
%               SOA_STEP_: 2740
%               FixError_: 2750
%              GoSaccade_: 2751
%                GoError_: 2752
%              NOGOWrong_: 2753
%         GoTargFixError_: 2754
%              GOCorrect_: 2755
%            NOGOCorrect_: 2756
%           CatchCorrect_: 2757
%        CatchIncorrectG_: 2758
%       CatchIncorrectNG_: 2759
%              BreakTFix_: 2760
%           EarlySaccade_: 2761
%              FixWindow_: 2770
%           TargetWindow_: 2771
%              Staircase_: 2772
%      Neg2Reinforcement_: 2773
%               Feedback_: 2774
%            ExtraReward_: 2777
%          SoundOnReward_: 2778
%          SoundNoReward_: 2779
%                Saccade_: 2810
%                 Decide_: 2811
%             RewardSize_: 2927
%            TrialInBlock: 2928
%         SendPenatrInfo_: 2929
%             StartInfos_: 2998
%               EndInfos_: 2999