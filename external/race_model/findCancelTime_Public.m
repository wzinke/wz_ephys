%% findCancelTime_Public.m

plotAT = 0;     % set to 1 to plot traces

%=============================================================
% on average traces - find finish time using average
% calculated from 1:SSD+Dstop
%=============================================================

numGreaterThanHI = 6;  % 6 times SD of noise
numGreaterThanLO = 2;  % 2 times SD of noise
timeGreaterThan = 50;   % amount of time numGreaterThanLO after numGreaterThanHI
for x = 1:1:length(SSDall)
    % canceled
    go = 1; f = 1;
    diffArr=[]; crossStd = [];
    diffArr = ActFunc(x).Can - ActFunc(x).NSSc;
    aveDiff = nanmean(diffArr(1:SSDall(x)));
    stdDiff = nanstd(diffArr(1:SSDall(x)));

    crossStd = find(abs(diffArr)>stdDiff*numGreaterThanLO);
    if ~isempty(crossStd)
        while (go ~= 0) & (f<length(crossStd))
            if (ActFunc(x).Can(crossStd(f)) > threshMove1) | (ActFunc(x).NSSc(crossStd(f)) > threshMove1)
                CancelTime_AT(1,x) = NaN;
                go = 0;
            else
                if (crossStd(f)+timeGreaterThan-1) < length(diffArr)
                    arrBig = []; arrP = [];
                    arrP = abs(diffArr(crossStd(f):crossStd(f)+timeGreaterThan-1));
                    arrBig = find(arrP>stdDiff*numGreaterThanHI);
                    if length(arrBig) == timeGreaterThan
                        CancelTime_AT(1,x) = crossStd(f);
                        go = 0;
                    else
                        f = f+1;
                        CancelTime_AT(1,x) = NaN;
                    end
                else
                    CancelTime_AT(1,x) = NaN;
                    go = 0;
                end
            end
        end
    else
        CancelTime_AT(1,x) = NaN;
    end

    % plot
    if plotAT == 1
        figure(88); subplot(2,2,x);
        plot(ActFunc(x).Can, 'r'); hold on;
        plot(ActFunc(x).NSSc, 'b'); hold on;
        plot(diffArr, 'k'); hold on;
        xlim([30 350])
        ylim([-200 1100])
        plot([SSDall(x) SSDall(x)], [0 1050], 'k', 'linewidth', 2); hold on;
        plot([CancelTime_AT(1,x) CancelTime_AT(1,x)], [0 1050], '--k', 'linewidth', 2); hold on;
        plot([0 350], [aveDiff aveDiff], 'k', 'linewidth', 2); hold on;
        plot([0 350], [aveDiff+stdDiff*numGreaterThanHI aveDiff+stdDiff*numGreaterThanHI], 'k:', 'linewidth', 2); hold on;
        plot([0 350], [aveDiff-stdDiff*numGreaterThanHI aveDiff-stdDiff*numGreaterThanHI], 'k:', 'linewidth', 2); hold on;
        if x == length(SSDall)
            legend('Can', 'NSS', 'diff', 'SSD', 'cancel time')
        end
    end

    % non-canceled
    go = 1; f = 1;
    diffArr=[]; crossStd = [];
    diffArr = ActFunc(x).NC - ActFunc(x).NSSnc;
    aveDiff = nanmean(diffArr(1:SSDall(x)));
    stdDiff = nanstd(diffArr(1:SSDall(x)));

    crossStd = find(abs(diffArr)>stdDiff*numGreaterThanLO);
    if ~isempty(crossStd)
        while (go ~= 0) & (f<length(crossStd))
            if (ActFunc(x).NC(crossStd(f)) > threshMove1) | (ActFunc(x).NSSnc(crossStd(f)) > threshMove1)
                CancelTime_AT(2,x) = NaN;
                go = 0;
            else
                if (crossStd(f)+timeGreaterThan-1) < length(diffArr)
                    arrBig = []; arrP = [];
                    arrP = abs(diffArr(crossStd(f):crossStd(f)+timeGreaterThan-1));
                    arrBig = find(arrP>stdDiff*numGreaterThanHI);
                    if length(arrBig) == timeGreaterThan
                        CancelTime_AT(2,x) = crossStd(f);
                        go = 0;
                    else
                        f = f+1;
                        CancelTime_AT(2,x) = NaN;
                    end
                else
                    CancelTime_AT(2,x) = NaN;
                    go = 0;
                end
            end % if passed thresh
        end % while loop
    else
        CancelTime_AT(2,x) = NaN;
    end % of isempty(crossStd)


    % plot
    if plotAT == 1
        figure(89); subplot(2,2,x);
        plot(ActFunc(x).NC, 'g'); hold on;
        plot(ActFunc(x).NSSnc, 'b'); hold on;
        plot(diffArr, 'k'); hold on;
        xlim([30 350])
        ylim([-200 1100])
        plot([SSDall(x) SSDall(x)], [0 1050], 'k', 'linewidth', 2); hold on;
        plot([CancelTime_AT(2,x) CancelTime_AT(2,x)], [0 1050], '--k', 'linewidth', 2); hold on;
        plot([0 350], [aveDiff aveDiff], 'k', 'linewidth', 2); hold on;
        plot([0 350], [aveDiff+stdDiff*numGreaterThanHI aveDiff+stdDiff*numGreaterThanHI], 'k:', 'linewidth', 2); hold on;
        plot([0 350], [aveDiff-stdDiff*numGreaterThanHI aveDiff-stdDiff*numGreaterThanHI], 'k:', 'linewidth', 2); hold on;
        if x == length(SSDall)
            legend('Can', 'NSS', 'diff', 'SSD', 'cancel time')
        end
    end % plot
end  % of SSD
% 'Finished average traces cancellation time estimate.    '


%==========================================================================
%           STOP PROCESS
%==========================================================================

%================================================
% COMPARE CAN AND NSS WITH A T-TEST
%================================================
STOPCancelTime_TT(1:length(SSDall)) = NaN;

for x=1:1:length(SSDall)

    % STOP process
    stopTtest = 0; i=1;timeCriterion = timeGreaterThan;
    while(stopTtest==0 & timeCriterion>0 & (SSDall(x)+time_delay+i)<=trialLength)

        CanArr = []; NSSCanArr = [];
        CanArr = ActFunc(x).STOPc(SSDall(x)+time_delay:SSDall(x)+time_delay+i);

        result = ttest_LB(CanArr);  % same as ttest, but nanmean function was different in .m file

        if result == 1
            if timeCriterion==timeGreaterThan
                STOPCancelTime_TT(x) = (SSDall(x)+time_delay+i);
            end
            timeCriterion = timeCriterion-1;
        else
            timeCriterion=timeGreaterThan;
        end

        if timeCriterion == 0
            stopTtest = 1;
        end

        i=i+1;
        if i>trialLength
            stopTtest = 1;
            STOPCancelTime_TT(x) = NaN;
        end
    end

 end

allCancelTimes.AT = ceil(CancelTime_AT);
allCancelTimes.TTstop = ceil(STOPCancelTime_TT);
