%% AccumulateResp_Public.m
% Simulate movement and fix cells (i.e., GO and STOP units)
%   activation(i) = [feedforward_input - decay + recurr_excitation - inhib_influences] + noise

% =========================================================================
%                               set TIME sections
% =========================================================================
endTime1 = fixRespond;                   % essentially only time 0
endTime2 = TrGO;                         % time = TrGO (accumulation in the GO process starts)
 
% stop signal presented
if aSTOPtrial == 1      
    endTime3 = SSDRespond+time_delay;    % time = SSD+Dstop (accumulation in the STOP process starts)
    endTime4 = trialLength;              % time = end of trial
else
    endTime3 = trialLength;              % no stop signal trial - accumulation in STOP never starts
    endTime4 = trialLength+1;
end

% figure out of what gets accumulated at each step (input to the GO and STOP units)
randDist1 = []; randDist1 = (normrnd(meanGO, stdGO, trialLength, 1));       % GO input
randDist2 = []; randDist2 = (normrnd(meanSTOP, stdSTOP, trialLength, 1));   % STOP input
stopwithGO = 0; stopwithSTOP = 0;

for trlen = 2:1:trialLength
        
    % pick a random increment for each time step (inter-trial variabilty)
    level_pGo = randDist1(trlen); 
    level_pStop = randDist2(trlen);
    
    % Determine if GO and/or STOP is active: if accGO/accStop set to 0, no
    % accumulation; set to 1 to accumulate
    if trlen == 1
    elseif trlen <= endTime1    % in time1 time period
        accGO = 0;
        accStop = 0;
    elseif trlen <= endTime2    % in time2 time period (before TrGO)
        accGO = 0;
        accStop = 0;
    elseif trlen <= endTime3    % in time3 time period (after TrGO, before SSD+Dstop)
        accGO = 1;
        accStop = 0;
    elseif trlen <= endTime4    % in time4 time period (after SSD+Dstop)
        accGO = 1;
        accStop = 1;
    end
    
    % current level of activity
    currentGo = move1cell(trlen-1);
    currentStop = fixcell(trlen-1);       
    
    % added noise 
    % NOTE: We calculated noise using vector approach (randDist1, randDist2), 
    % instead of trial-by-trial approach (setting noiseGo and noiseStop),
    % because this was a bit quicker and produced the same results
    noiseGo = 0;
    noiseStop = 0;
    
    % GO unit activation:
    move1cell(trlen) = ((currentGo*accGO + level_pGo*accGO) - ...
        (kGo*currentGo) - ((BStop*accStop)*currentStop)) + noiseGo;
    % STOP unit activation:
    fixcell(trlen) = ((currentStop*accStop + level_pStop*accStop) - ...
        (kStop*currentStop) - ((BGo*accGO)*currentGo)) + noiseStop;

    % rectify activations if they fall below 0    
    if move1cell(trlen) < 0
        move1cell(trlen) = 0;
    end
    if fixcell(trlen) < 0 
        fixcell(trlen) = 0;
    end 

end