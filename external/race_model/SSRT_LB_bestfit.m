function [meanIntSSRT, meanSSRT, overallMeanSSRT, HanesSSRT, meanInhibFunc, Time, BestWFit, FitValues, meanIntSSRT_bf] = ...
    SSRT_LB_bestfit(NSSDist, inhib_func, SSDall)
%% SSRT_LB_bestfit.m

% written by L. Boucher March 2003
%======================================================================
% INPUTS:
% NSSDist = array with all no stop signal RTs
% inhib_func = inhibition function
% SSDall = array with all SSDs
%======================================================================
% OUTPUTS:
% meanIntSSRT = via integration method
% meanSSRT = via difference method
% overallMeanSSRT = average of both
% Time = time ticks from lowest SSD to highest
% BestWFit = best fit weibell value of inhibition function
% FitValues = fit values from weibell mapped from SSDs, x value for BestWFit values
% meanIntSSRT_bf = mean of the int method when calculated from weibell
%======================================================================

%======================================================================
% rank NSS distribution
%======================================================================
if ~isempty(NSSDist)
    meanNSS = nanmean(NSSDist);
    stdNSS = nanstd(NSSDist);
    NNSS = length(find(~isnan(NSSDist)));
else
    meanNSS = NaN;
    stdNSS = NaN;
    NNSS = NaN;
end

% rank and get rid of NaNs in RT dist
rankedNSS_RT = sort(NSSDist);
rankedNSS_RT = rankedNSS_RT(find(~isnan(rankedNSS_RT)));

%======================================================================
% calculate mean SSRT after Logan and Cowan (1984)
% difference between mean RT and mean oif inhib func (as calculated 
% using the Weibell best fit function). 
%======================================================================
% get SSRT by Weibull (assuming it to be a random variable)
% get initial estimate
[d,ClosestSSD] = min(abs(inhib_func-0.67));
InitialGuess(1) = SSDall(ClosestSSD); % time inhib func reaches 67% of its max
InitialGuess(2) = 10;                 % slope
InitialGuess(3) = 1;                  % max of inhib func 
InitialGuess(4) = 0;                  % min of inhib func

% get Non-Linear Fit using Least-Square Fit

Params = nlinfit(SSDall,inhib_func,'InhCumWeib',InitialGuess);

Time = [round(min(SSDall)):1:round(max(SSDall))];

BestWFit = InhCumWeib(Params,Time);

FitValues = BestWFit(SSDall-(min(SSDall))+1);

for x = 2:1:(length(Time)-1)
    val(x-1) = (BestWFit(x) - BestWFit(x-1)) * Time(x);
end
meanInhibFunc = sum(val)/(max(FitValues) - min(FitValues));
meanSSRT = meanNSS - meanInhibFunc;


%======================================================================
% calculate SSRTs after Hanes et al (1998) - INTEGRATION method (originally Logan
% and Cowan, 1984)
% 1. rank order NSS RTs
% 2. mulitply the p(noncan) by total number of NSS trials to get
%    index 
% 3. find ranked ordered NSS at index above and subtract SSD
% 
% USE INHIBITION FUNCTION 
%======================================================================

for x = 1:1:length(SSDall)
    indSSRT(x) = round(inhib_func(x)*NNSS);
    if (isnan(indSSRT(x)) == 0) & (indSSRT(x) ~=0)
        HanesSSRT(x) = rankedNSS_RT(indSSRT(x)) - SSDall(x);
    else
        HanesSSRT(x) = NaN;
    end
end
meanIntSSRT = nanmean(HanesSSRT);


%======================================================================
% calculate SSRTs after Hanes et al (1998) - INTEGRATION method
% 1. rank order NSS RTs
% 2. mulitply the p(noncan) by total number of NSS trials to get
%    index 
% 3. find ranked ordered NSS at index above and subtract SSD
%
% USE BEST FIT VALUES
%======================================================================
FitValues_bounded = FitValues;
ind = [];
ind = find(FitValues>1);
FitValues_bounded(ind) = 1;
ind = [];
ind = find(FitValues<0);
FitValues_bounded(ind) = 0;

for x = 1:1:length(SSDall)
    indSSRT_bf(x) = round(FitValues_bounded(x)*NNSS);
    if isreal(indSSRT_bf) ~= 0
        if (isnan(indSSRT_bf(x)) == 0) & (indSSRT_bf(x) ~=0)
            HanesSSRT_bf(x) = rankedNSS_RT(indSSRT_bf(x)) - SSDall(x);
        else
            HanesSSRT_bf(x) = NaN;
        end
    else
        HanesSSRT_bf(x) = NaN;
    end
end
meanIntSSRT_bf = nanmean(HanesSSRT_bf);


%======================================================================
% Overall Stats
%======================================================================
overallMeanSSRT = nanmean([meanSSRT meanIntSSRT]);

% round to 2 decimal places
meanIntSSRT = (ceil(meanIntSSRT*1000))/1000;
meanSSRT = (ceil(meanSSRT*1000))/1000;
overallMeanSSRT = (ceil(overallMeanSSRT*1000))/1000;
