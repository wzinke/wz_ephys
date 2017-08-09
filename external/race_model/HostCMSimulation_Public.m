function [NSS_rtSim, inhib_funcSim, SS_rtSim, predSSRT, allCancelTimes, ActFunc] = ...
    HostCMSimulation_Public(param, SSDall, NUMB_TR, modelNo, master_seed, iteration);
%% HostCMSimulations_Public.m

warning off MATLAB:divideByZero
threshMove1 = 1000;     % for GO process
threshFix = 1000;       % for STOP process
trialLength = 1200;     % time simulation works over
numObs = NUMB_TR;
% because the NSS are compared with both canceled and noncanceled trials,
% more NSS trials need to be run in the sub-sample iteration for
% calculation of cancel time
numObs_NSS = (numObs*length(SSDall))+(numObs*length(SSDall)*2);
maxNumObs = ceil(numObs*2.5);   % will only simulate a maximum of this many trials
numInBin = 0; iSSD = 1; addme = 0;
TOTnumInBin = 0; STOPnow = 0; new = 0;

if iteration == 0
    numTrials = numObs + maxNumObs*length(SSDall);      % (numObs for NSS trials; numObs for each SSD)
else
    numTrials = numObs_NSS + maxNumObs*length(SSDall);  % (numObs for NSS trials; numObs for each SSD)
end

%==========================================================
% initialize arrays
%==========================================================
goProcess(1:numTrials, 1:trialLength) = NaN;
stopProcess(1:numTrials, 1:trialLength) = NaN;
nanarr(1:trialLength) = nan;

%==========================================================
% run a bunch of trials
%==========================================================
nt = 0; nnt = 0;maxNumTrials=numTrials; t = 1;
while (nt < maxNumTrials & STOPnow == 0)
    nt = nt+1;
    randn('state',nt+master_seed);
    rand('state',nt+master_seed);

    %==========================================================
    % set parameters
    %==========================================================
    SetParams_Public;

    %==========================================================
    % run one trial of simulation
    %==========================================================
    OneTrial_Public;

    STOPtrial(nt) = aSTOPtrial;
    move1RT(nt) = amove1RT;
    fixRT(nt) = afixRT;
    SSD(nt) = aSSD;

    goProcess(nt, :)  = move1cell;
    stopProcess(nt, :) = fixcell;

    % =====================================================================
    % DECISION RACE (independent model):
    % 1 = go; 0 = stop
    % if process1 doesn't have an RT, and process2 does, process2 wins
    % if neither process has an RT, decision is to STOP
    % =====================================================================
    if (isnan(move1RT(nt)) == 0) & (isnan(fixRT(nt)) == 0)  % neither equal to NaN
        if move1RT(nt) < fixRT(nt);
            decRACE(nt) = 1;   % GO
        else move1RT(nt) > fixRT(nt);
            decRACE(nt) = 0;   % STOP
        end
    elseif (isnan(move1RT(nt)) == 1) & (isnan(fixRT(nt)) == 0); % moveRT equal to NaN
        decRACE(nt) = 0;   % STOP
    elseif (isnan(move1RT(nt)) == 0) & (isnan(fixRT(nt)) == 1); % fixRT equal to NaN
        decRACE(nt) = 1;   % GO
    elseif (isnan(move1RT(nt)) == 1) & (isnan(fixRT(nt)) == 1); % both equal to NaN
        decRACE(nt) = 0;   % STOP
    end

    % =====================================================================
    % DECISION GO (interactive model):
    % 1 = go hits threshold; 0 = go DOES NOT hit threshold
    % =====================================================================
    if ~isnan(move1RT(nt))
        decGO(nt) = 1;
    else
        decGO(nt) = 0;
    end

    % =====================================================================
    % DECISION is based on which rule?
    % 1 = decGO, non-canceled is if GO hits threshold
    % 0 = decRACE, nan-canceled is if RTgo < RTstop
    % =====================================================================
    if modelNo == 1
        decision(nt) = decGO(nt);       % interactive rule
    else
        decision(nt) = decRACE(nt);     % independent rule
    end

    % =====================================================================
    % ADD 10ms to RTgo to account for time needed to move the eyes
    % =====================================================================
    ballisticTime = 10;
    move1RT(nt) = move1RT(nt) + ballisticTime;

end

%==========================================================
% define NSS trials:
%==========================================================
NSS_ind = find(STOPtrial==0);
NSS_rtSim = move1RT(NSS_ind);

%==========================================================
% Calculate inhibition function and SSRT for each SSD
% %==========================================================
CalcSSRT_Public;

%==========================================================
% calculate SSRT
%==========================================================
predSSRT(1)= meanIntSSRT;
predSSRT(2) = meanSSRT;
predSSRT(3) = overallMeanSSRT;
predSSRT(4:3+length(SSDall)) = SSRT_by_SSD;

%==========================================================
% find latency matched trials
%==========================================================
LatMatch_Public;

%==========================================================
% get model activation functions for canceled, non-canceled, NSS trials, and STOP process
%==========================================================
getActFunc_Public;


%==========================================================
% find cancel times for GO and STOP unit.
%==========================================================
findCancelTime_Public;
