%% LatMatch_Public.m
ActFunc = []; CancelTime=[];

Divider = targON + SSDall + overallMeanSSRT;

% ****************************************************
% 1. finds indices above and below SSRT
% ****************************************************

for ev = 1:1:length(SSDall)+1
    if ev == 1
        for x = 1:1:length(SSDall)
            dd = []; dd = find(NSS_rtSim > Divider(x));
            NSS_lat(x).above = NSS_ind(dd);
            dd = []; dd = find(NSS_rtSim <= Divider(x));
            NSS_lat(x).below = NSS_ind(dd);
        end
    else
        dd = []; dd = find(SS_rtSim(ev-1).dec == 0);    % canceled
        SS_lat(ev-1).above = SS_rtSim(ev-1).ind(dd);
        dd = []; dd = find(SS_rtSim(ev-1).dec == 1);    % non-canceled
        SS_lat(ev-1).below = SS_rtSim(ev-1).ind(dd);
    end
end
