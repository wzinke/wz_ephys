function [Tmat, plxMat] = plx_trialtable(plxMat, tblfl, EV, NHP, ITIstim, savemat)

%EV ='/DATA/wz_code/plexon/TEMPO_EV_cosman_rig028.m';
rig=28;

% ____________________________________________________________________________ %
%% check input data
% if not specified use GUI to get the file
if(~exist('plxMat','var') || isempty(plxMat))
    [FileName,PathName] = uigetfile({'*.plx;*.PLX'},'Load plexon file');
    plxMat = fullfile(PathName,FileName);
end

if(~exist('ITIstim','var') || isempty(ITIstim))
    ITIstim = 0; % flag to indicate whether stimulation occurred between trial (ad hoc hack for ultrasound experiments)
end


if(~isstruct(plxMat))
    % get event code information
    if(~exist('EV','var') || isempty(EV))
        [FileName,PathName] = uigetfile({'*.pro;*.mat;*.m'},'Event code definition file');
        EV = fullfile(PathName,FileName);
    end

    if(~exist('NHP','var') || isempty(NHP))
        if(isstruct(plxMat))
            NHP = plxMat.File(1);
        else
            [~,cnm] = fileparts(plxMat);
            NHP = cnm(1); % assume default that first letter specifies NHP
        end
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

    %% read data file
    plxMat = PLX_get_paradigm(plxMat, EV, rig, NHP, 0, 1, 1, ITIstim);
end

% standard info to write
trialtime = plxMat.TrialStartEndDur(:,1) + plxMat.Task.StimOnsetToTrial;
vp = isfinite(trialtime);

TBL.NHP       = repmat(plxMat.NHP,sum(vp),1);
TBL.Date      = repmat(datestr(plxMat.DateNum,'yyyy_mm_dd'),sum(vp),1);
TBL.TrialNum  = [1:sum(vp)]';
TBL.TrialTime = plxMat.TrialStartEndDur(vp,1) + plxMat.Task.StimOnsetToTrial(vp);
TBL.Task      = plxMat.Task.TaskType(vp)';
TBL.Correct   = plxMat.Task.Correct(vp);
TBL.Error     = plxMat.Task.error(vp);
TBL.SRT       = plxMat.Task.SRT(vp);

% list of extra info
fldlst = sort({'TargetHemi', 'TargetPos', 'TargetLoc', 'SetSize', 'DistHomo', 'SingMode', 'TargetType', ...
          'Eccentricity', 'Singleton', 'DistHemi', 'DistLoc', 'DistPos', 'IsCatch', 'StimTrial', ...
          'DelayDur', 'MGcolor', 'SingleDistCol', 'TargetOri', 'DistOri', 'FixSpotOff', 'FixSpotOn', ...
          'SOA', 'StimTime','StimOnsetToTrial', 'StimTrialTime', 'StimBlockTrial', 'StimBlockCond', ...
          'SaccStart', 'SaccEnd', 'SaccDur', 'SacAng', 'SaccAmp', 'SaccSpeed', 'SaccDist', ...
          'Startpoint_X', 'Startpoint_Y', 'Endpoint_X', 'Endpoint_Y', ...
          'fix_X', 'fix_Y', 'fix_CI_Xrng', 'fix_CI_Yrng', ...
          'pre_X', 'pre_Y', 'pre_CI_Xrng', 'pre_CI_Yrng', ...
          'post_X post', 'Y post_CI_Xrng', 'post_CI_Yrng'});

for(i=1:length(fldlst))
    cfld = cell2mat(fldlst(i));

    if(isfield(plxMat.Task,cfld))
%        if(length(unique(plxMat.Task.(cfld))) > 1)
           TBL.(cfld) = plxMat.Task.(cfld)(vp);
%        end
    end
end

% convert to table and write it
Tmat = struct2table(TBL);
writetable(Tmat,tblfl);

if(exist('savemat','var') && ~isempty(savemat) && savemat == 1)
    [flpth, flstem] = fileparts(tblfl);
    save([flpth, '/', flstem, '.mat'], '-struct', 'plxMat', '-v7.3');
end



