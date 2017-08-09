function plx = PLX_EyeTrace(plx)
% Use a data driven approach to determine SRT
%
% DESCRIPTION
%       
%
% SYNTAX
%
%   plx = PLX_EyeTrace(plx)
%
%   Input:
%
%
%
% .........................................................................
% wolf zinke, wolfzinke@gmail.com
%
% $Created : 12-Jul-2015 by wolf zinke

% initialize new structure fields

plx.Eye.smX    = nan(size(plx.Eye.EyeX));
plx.Eye.smY    = nan(size(plx.Eye.EyeY));
plx.Eye.smDist = nan(size(plx.Eye.EyeX));

plx.Task.SaccStart = nan(plx.NTrials,1);
plx.Task.SaccEnd   = nan(plx.NTrials,1);
plx.Task.SaccDur   = nan(plx.NTrials,1);
plx.Task.SaccAng   = nan(plx.NTrials,1);
plx.Task.SaccAmp   = nan(plx.NTrials,1);
plx.Task.SaccSpeed = nan(plx.NTrials,1);
plx.Task.SaccDist  = nan(plx.NTrials,1);

plx.Task.Startpoint_X = nan(plx.NTrials,1);
plx.Task.Startpoint_Y = nan(plx.NTrials,1);
plx.Task.Endpoint_X   = nan(plx.NTrials,1);
plx.Task.Endpoint_Y   = nan(plx.NTrials,1);

plx.Task.fix_X = nan(plx.NTrials,1);
plx.Task.fix_Y = nan(plx.NTrials,1);
plx.Task.fix_CI_Xrng = nan(plx.NTrials,1);
plx.Task.fix_CI_Yrng = nan(plx.NTrials,1);

plx.Task.pre_X = nan(plx.NTrials,1);
plx.Task.pre_Y = nan(plx.NTrials,1);
plx.Task.pre_CI_Xrng = nan(plx.NTrials,1);
plx.Task.pre_CI_Yrng = nan(plx.NTrials,1);

plx.Task.post_X = nan(plx.NTrials,1);
plx.Task.post_Y = nan(plx.NTrials,1);
plx.Task.post_CI_Xrng = nan(plx.NTrials,1);
plx.Task.post_CI_Yrng = nan(plx.NTrials,1);

for(i=1:plx.NTrials)
    if(isnan(plx.Task.StimOnset))
        continue;
    end
    % disp(i)
    Gaze = SRT_EyeTraces(plx.Eye.EyeX(i,:), plx.Eye.EyeY(i,:), plx.timevec, plx.Task.StimOnset(i), plx.Task.GoCue(i), 0);

    plx.Eye.smX(i,:)    = Gaze.smX;
    plx.Eye.smY(i,:)    = Gaze.smY;
    plx.Eye.smDist(i,:) = Gaze.smD;

    plx.Task.SaccStart(i) = Gaze.SaccStart;
    plx.Task.SaccEnd(i)   = Gaze.SaccEnd;
    plx.Task.SaccDur(i)   = Gaze.SaccDur;
    plx.Task.SaccAng(i)   = Gaze.SaccAng;
    plx.Task.SaccAmp(i)   = Gaze.SaccAmp;
    plx.Task.SaccSpeed(i) = Gaze.SaccSpeed;
    plx.Task.SaccDist(i)  = Gaze.SaccDist;

    plx.Task.Startpoint_X(i) = Gaze.Startpoint_X;
    plx.Task.Startpoint_Y(i) = Gaze.Startpoint_Y;
    plx.Task.Endpoint_X(i)   = Gaze.Endpoint_X;
    plx.Task.Endpoint_Y(i)   = Gaze.Endpoint_Y;

    plx.Task.fix_X(i) = Gaze.fix_X;
    plx.Task.fix_Y(i) = Gaze.fix_Y;
    plx.Task.fix_CI_Xrng(i) = Gaze.fix_CI_Xrng;
    plx.Task.fix_CI_Yrng(i) = Gaze.fix_CI_Yrng;

    plx.Task.pre_X(i) = Gaze.pre_X;
    plx.Task.pre_Y(i) = Gaze.pre_Y;
    plx.Task.pre_CI_Xrng(i) = Gaze.pre_CI_Xrng;
    plx.Task.pre_CI_Yrng(i) = Gaze.pre_CI_Yrng;

    plx.Task.post_X(i) = Gaze.post_X;
    plx.Task.post_Y(i) = Gaze.post_Y;
    plx.Task.post_CI_Xrng(i) = Gaze.post_CI_Xrng;
    plx.Task.post_CI_Yrng(i) = Gaze.post_CI_Yrng;
end