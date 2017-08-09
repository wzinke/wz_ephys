function [TDT001 , TDT005] = wz_spk_respdiff_fast(spk1, spk2, timwin, nBoot, mincons, siglev)
% Determines the time where two response profiles start to become different
% from each other. 
%
% wolf zinke, 22.1.2014


%% check input and define default values
if(~isfield(spk1,'spikedensities'))
    spk1 = wz_spk_density(spk1, 'gauss', 20);
end

if(~isfield(spk2,'spikedensities'))
    spk2 = wz_spk_density(spk2, 'gauss', 20);
end

% ensure a time window is used for the test that suits both series 
if(isempty(spk1.timewindow) || isempty(spk2.timewindow))
    warning('At least one of the spike matrices do not contain any data! \n ... nothing to compare here.');
    TDT005 = NaN;
    TDT001 = NaN;
    return
end

if(~exist('timwin','var') || isempty(timwin))
    timwin(1) = max([spk1.timewindow(1); spk2.timewindow(1)]);
    timwin(2) = min([spk1.timewindow(2); spk2.timewindow(2)]);
elseif(any([spk1.timewindow(1);spk2.timewindow(1)] > timwin(1)) || ...
   any([spk1.timewindow(2);spk2.timewindow(2)] < timwin(2)) )
    timwin(1) = max([spk1.timewindow(1); spk2.timewindow(1)]);
    timwin(2) = min([spk1.timewindow(2); spk2.timewindow(2)]);
end

% minimum requirement of consecutive significant time bins
if(~exist('nBoot','var') || isempty(nBoot))
    nBoot = 1;
end

% minimum requirement of consecutive significant time bins
if(~exist('mincons','var') || isempty(mincons))
    mincons = 10;
end


tvec = timwin(1):timwin(2);  % time bins used for the test

%% get the trial based spike densities
tp1 = spk1.xtime >= timwin(1) & spk1.xtime <= timwin(2);
tp2 = spk2.xtime >= timwin(1) & spk2.xtime <= timwin(2);

spkmat1 = spk1.spikedensities(:, tp1);
spkmat2 = spk2.spikedensities(:, tp2);

% get the boot sample matrix
rng(sum(100*clock),'twister');

Bsmpl1 = randi(spk1.nTrials, spk1.nTrials, nBoot);
Bsmpl1(:,1) = 1:spk1.nTrials; % make sure that first sample is original data
Bsmpl2 = randi(spk2.nTrials, spk2.nTrials, nBoot);
Bsmpl2(:,1) = 1:spk2.nTrials;

%% run bootstrap loops 
TDT001 = nan(1,nBoot);
TDT005 = nan(1,nBoot);

for(r=1:nBoot)
    % get the current trial sample (with replacement)   
    cspikes1 = spkmat1(Bsmpl1(:,r),:);
    cspikes2 = spkmat2(Bsmpl2(:,r),:);
    
    pvals    = nan(1, length(tvec));
    
    sigcnt = 0;
    for(b=1:length(tvec))
        csmpl1 = cspikes1(:,b);
        csmpl2 = cspikes2(:,b);
%         csmpl1(isnan(csmpl1)) = [];
%         csmpl2(isnan(csmpl2)) = [];

%         if(~isempty(csmpl1) || ~isempty(csmpl2))
            % use bin wise rank sum test to find response differences
            [pvals(b)] = ranksum(csmpl1, csmpl2, 'method','approximate');
%         else
%             sigcnt = 1; 
%         end
        
        if(pvals(b) < 0.01)
            sigcnt = sigcnt + 1;
        else
            sigcnt = 0;
        end
        
        if(sigcnt == mincons) % citerion fullfilled
            ct = b-mincons+1;
            TDT001(r) = tvec(ct);
            
            pos = find(pvals(1:ct) < 0.05,1,'last');
            if(~isempty(pos))
                TDT005(r) = tvec(pos)+1;            
            end
            
            break;
        end
    end
end

