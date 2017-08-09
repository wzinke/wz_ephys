%% CalcSSRT_Public.m
numSSD = length(SSDall);

for x = 1:1:numSSD+1
    ind = []; dec = []; stop = []; go = []; SSRT = [];
    %first time around, get NSS data sorted out
    if x==numSSD+1
        ind = find(SSD==0);
        NSS_stop = fixRT(ind);
        simNSS_rt = move1RT(ind);
        NSS_ind = ind;
    else
        ind = find(SSD==SSDall(x));

        dec = decision(ind);
        stop = fixRT(ind);
        go = move1RT(ind);

        SS_rtSim(x).rt = move1RT(ind);
        SS_rtSim(x).ind = ind;
        SS_rtSim(x).dec = dec;

        prop_noncan = sum(dec)/(length(dec));
        inhib_funcSim(x) = prop_noncan;

        SSRT = mean(stop)-SSDall(x);
        SSRTmean(x) = SSRT;
    end
end

[meanIntSSRT, meanSSRT, overallMeanSSRT, SSRT_by_SSD ]= ...
    SSRT_LB_bestfit(simNSS_rt, inhib_funcSim, SSDall);
