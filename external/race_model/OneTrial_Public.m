%% OneTrial_Public.m
targON = 0; % target is presented at time 0

%=========================================================================
% determine # in each bin.  Keep going until a minimum number of RTs is
% reached (numObs) or a maximum of trials is run in that case (maxNumObs)
%=========================================================================
numInBin = numInBin + addme;

if iteration == 0
    NO = numObs;
else
    NO = numObs_NSS;
end
    
if nt <= NO
    aSTOPtrial = 0;  % no STOP signal
    aSSD = 0;
else
    aSTOPtrial = 1;  % STOP signal
    
    if TOTnumInBin < (maxNumObs-1) & numInBin < (numObs-1) 
        aSSD = SSDall(iSSD);
        TOTnumInBin = TOTnumInBin+1;
        new = 1;
    else
        new = 0;
        aSSD = SSDall(iSSD);
        iSSD = iSSD+1;
        if iSSD == length(SSDall)+1
            STOPnow = 1;
        end
        numInBin = 0;
        TOTnumInBin = 0;
    end
end

%===========================================================
% time fixation takes to get into the system relative to trial start
% (this is actually TrGO)
%===========================================================
fixRespond = 0;
fixActual = 0;

%===========================================================
% time target gets into the system relative to trial start
%===========================================================
if targON == 0
    targONa = 1;
end

targRespond = (targONa);
targActual = targON;

%===========================================================
% time Stop Signal gets into the system relative to trial start
%===========================================================
if aSTOPtrial == 1
    SSDRespond = (targON + aSSD);
    SSDActual = (targON + aSSD);
else
    SSDRespond = 0;
    SSDActual = 0;
end

% initialize arrays with NaNs
move1cell = []; move1cell(1:trialLength) = NaN; move1cell(1) = 0;   % GO unit
fixcell = []; fixcell(1:trialLength) = NaN; fixcell(1) = 0;         % STOP unit

% accumulate activation in the GO and STOP units
AccumulateResp_Public; 
       
% GO and STOP RT measured relative to beginning of trial, after targ onset
moveth = []; fixth = [];
moveth = (find(move1cell >= threshMove1));
fixth = (find(fixcell >= threshFix));

% find RT on trial relative to targ onset
if length(moveth) >=1
    amove1RT = moveth(1);
else
    amove1RT = NaN;
end
if length(fixth) >=1
    afixRT = fixth(1);
else
    afixRT = NaN;
end

% for bookkeeping number of trials
if aSTOPtrial==1 & new ==1
    if decision(nt-1)==1
        addme=1;
    else
        addme=0;
    end
end

