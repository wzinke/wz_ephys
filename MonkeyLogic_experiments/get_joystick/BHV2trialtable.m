function [Otbl, obhv] = BHV2trialtable(ifl, ofl)
% This function translates a *.bhv file (<ifl>) created by MonkeyLogic into a simple trial table.
% This table could be saved as text file (<ofl>) to be further processed with R, or a struct as
% output (<Otbl>) could be used for further analysis within Matlab.
%
% This routine relies on definitions of event codes as they are used in the
% get_joystick.m timing file. Use the output generated with this function to
% check for consistency with the textfile that is created in this timing file.
%
% wolf zinke, Apr. 2014

%% make sure MonkeyLogic directory is known
if(exist('bhv_read','file') ~= 2)
    MLdir = uigetdir(pwd,'MonkeyLogic Directory');
    addpath(MLdir);
end

%% get data file
if(exist('ifl','var') == 0 || isempty(ifl) == 1)
    [ifl, pathname] = uigetfile([pwd '*.bhv'], 'Choose BHV file');
    if ifl == 0,
        return
    end
    ifl = [pathname ifl];
end

if(isstruct(ifl) == 1)
    bhv = ifl;
else
    bhv = bhv_read(ifl);
end

%% read in the data file
tpos      = find(bhv.TrialError ~= 9 & bhv.TrialError ~= 8); % ignore skipped trials
numTrials = length(tpos);

%% pre-define trialtable
Ttbl.Date = datestr(bhv.AbsoluteTrialStartTime(tpos,:),'YYYY_MM_DD');
Ttbl.StartTime = datestr(bhv.AbsoluteTrialStartTime(tpos,:),'hh:mm:ss');
Ttbl.Subject   = repmat(bhv.SubjectName,numTrials,1);

Ttbl.trialnumber = bhv.TrialNumber(tpos);
Ttbl.Cond        = bhv.ConditionNumber(tpos);

Ttbl.TrialError  = bhv.TrialError(tpos);
Ttbl.correct     = bhv.TrialError(tpos) == 0;
Ttbl.RT          = bhv.ReactionTime(tpos)';

for(t=1:numTrials)
    EvCodes = cell2mat(bhv.CodeNumbers(tpos(t)));
    EvTimes = cell2mat(bhv.CodeTimes(tpos(t)));
    
    %% translate event codes into meaningful names for variables
    FixOn      = Get_time(EvTimes,EvCodes, 35);
    FixOff     = Get_time(EvTimes,EvCodes, 36);
    TrialStart = Get_time(EvTimes,EvCodes,  1);
    TrialEnd   = Get_time(EvTimes,EvCodes,  2);
    DimmTime   = Get_time(EvTimes,EvCodes, 37); 
    RelTime    = Get_time(EvTimes,EvCodes, 39);
    RewEv      = min(Get_time(EvTimes,EvCodes, 96));
    RewTime    = bhv.RewardRecord(t).RewardOnTime;
    RewTime(RewTime<TrialEnd) = [];
    numRew = length(RewTime);
    
    %% add details to the trial table   
    Ttbl.WaitStart(t,1) = TrialStart - FixOn;
    Ttbl.StimDur(t,1)   = FixOff     - FixOn;
    Ttbl.DimmTime(t,1)  = DimmTime   - TrialStart;   
    Ttbl.RelTime(t,1)   = RelTime    - TrialStart;    
    if(Ttbl.TrialError(t) == 0)
        Ttbl.RTcalc(t,1) = RelTime    - DimmTime; % this should be the same as bhv.ReactionTime(tpos(t))!
    else
        Ttbl.RTcalc(t,1) = NaN;
    end
    Ttbl.WaitRew(t,1)   = RewEv      - RelTime;     
    
    Ttbl.numRew(t,1)    = numRew; 
end

%% create output
if(nargin == 2)
    TBL = struct2table(Ttbl);
    writetable(TBL,ofl);
end

if(nargout > 0)
   Otbl = Ttbl; 
end
if(nargout > 1)
   obhv = bhv; 
end

%% just simplify checking for event times
function ctm = Get_time(EvT, EvC, cEv)
    ctm = EvT(EvC == cEv);
    if(isempty(ctm))
        ctm = NaN;
    end

    
%% Event codes used in this paradigm
% strcat(num2str(bhv.CodeNumbersUsed), '----', bhv.CodeNamesUsed)
%     ' 1----Start trial'
%     ' 2----End trial'
%     ' 9----Start inter-trial interval'
%     '10----End inter-trial interval'
%     '11----Start wait fixation'
%     '12----End wait fixation'
%     '15----Start pre-trial'
%     '16----End pre-trial'
%     '17----Start post-trial'
%     '18----End post-trial'
%     '19----Start pause'
%     '20----End pause'
%     '35----Fixation spot ON'
%     '36----Fixation spot OFF'
%     '37----Dim fixation spot'
%     '38----Start up lever'
%     '39----End up lever'
%     '96----Reward delivered'



