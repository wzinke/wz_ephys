%% getActFunc_Public.m

% close all;
ActFunc = [];  boundActFunc=[];

NaNarray(1:trialLength) = NaN;

plotAllTrials = 0;
plotFill = 0;
plotSTOP = 1;

for x = 1:1:length(SSDall)
    %==========================================
    % CANCELLED
    %==========================================
    
    Xmean1 = [];Xcumsum=[];lowInd = []; highInd = [];
    NumActFunc(1, x) = length(SS_lat(x).above);  % Can, NSSc, NC, NSSnc
    if length(SS_lat(x).above)>1
        Xmean1 = nanmean(goProcess(SS_lat(x).above, :));
    else
        Xmean1 = (goProcess(SS_lat(x).above, :));
    end
    Xcumsum = sort(goProcess(SS_lat(x).above, :));
    if ~isempty(Xcumsum)
        ActFunc(x).Can = Xmean1;
    else
        ActFunc(x).Can = NaNarray;
    end
    
    
    Xmean2 = [];Xcumsum=[];lowInd = []; highInd = [];
        NumActFunc(2, x) = length(NSS_lat(x).above);  % Can, NSSc, NC, NSSnc
    if length(NSS_lat(x).above)>1
        Xmean2 = nanmean(goProcess(NSS_lat(x).above, :));
    else
        Xmean2 = (goProcess(NSS_lat(x).above, :));
    end
    Xcumsum = sort(goProcess(NSS_lat(x).above, :));
    if ~isempty(Xcumsum)
        ActFunc(x).NSSc = Xmean2;
    else
        ActFunc(x).NSSc = NaNarray;
    end
    
    if plotSTOP == 1
        Stopmean = [];Xcumsum=[];lowInd = []; highInd = [];
        if length(SS_lat(x).above)>1
            Stopmean = nanmean(stopProcess(SS_lat(x).above, :));
        else
            Stopmean = (stopProcess(SS_lat(x).above, :));
        end
        Xcumsum = sort(stopProcess(SS_lat(x).above, :));
        if ~isempty(Xcumsum)
            ActFunc(x).STOPc = Stopmean;
        else
            ActFunc(x).STOPc = NaNarray;
        end
        
    end
    
    %==========================================
    % NON-CANCELLED
    %==========================================
    Xmean1 = [];Xcumsum=[];lowInd = []; highInd = [];
            NumActFunc(3, x) = length(SS_lat(x).below);  % Can, NSSc, NC, NSSnc
    if length(SS_lat(x).below)>1
        Xmean1 = nanmean(goProcess(SS_lat(x).below, :));
    else
        Xmean1 = (goProcess(SS_lat(x).below, :));
    end
    Xcumsum = sort(goProcess(SS_lat(x).below, :));
    if ~isempty(Xcumsum)
        ActFunc(x).NC = Xmean1;
    else
        ActFunc(x).NC = NaNarray;
    end
    
    Xmean2 = [];Xcumsum=[];lowInd = []; highInd = [];
    NumActFunc(4, x) = length(NSS_lat(x).below);  % Can, NSSc, NC, NSSnc
    if length(NSS_lat(x).below)>1
        Xmean2 = nanmean(goProcess(NSS_lat(x).below, :));
    else
        Xmean2 = (goProcess(NSS_lat(x).below, :));
    end
    Xcumsum = sort(goProcess(NSS_lat(x).below, :));
    if ~isempty(Xcumsum)
        ActFunc(x).NSSnc = Xmean2;
    else
        ActFunc(x).NSSnc = NaNarray;
    end
    
    if plotSTOP == 1
        Stopmean = [];Xcumsum=[];lowInd = []; highInd = [];
        if length(SS_lat(x).below)>1
            Stopmean = nanmean(stopProcess(SS_lat(x).below, :));
        else
            Stopmean = (stopProcess(SS_lat(x).below, :));
        end
        
        Xcumsum = sort(stopProcess(SS_lat(x).below, :));
        if ~isempty(Xcumsum)
            ActFunc(x).STOPnc = Stopmean;
        else
            ActFunc(x).STOPnc = NaNarray;
        end
    end
    
end