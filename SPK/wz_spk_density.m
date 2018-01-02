function spk = wz_spk_density(spk, ktype, kwdth, rep, clip)
% Convolve a binary spike train with a kernel function in order to to
% obtain a response estimate.
% If a binary 2D matrix is used as input argument, the rows will be
% interpreted as trials. It is assumed that values represent spike times.
%
% wolf zinke, 13.1.2014
%
% ToDo: - Currently there seems to be an edge effect that results in flanking
%         response peaks. This needs to be debugged!

%% define default values
if(~exist('ktype','var') || isempty(ktype))
    ktype = 'gauss';
end

if(~exist('kwdth','var') || isempty(kwdth))
    kwdth = 10;
end

if(~exist('rep','var')   || isempty(rep))
    rep = 1;
end

if(~exist('clip','var')  || isempty(clip))
    doclip = 0;
else
    doclip = 1;
end

%% check data input
if(~isstruct(spk))
%     if(length(unique(spk(~isnan(spk)))) > 2)
        spk = wz_get_SPKobj(spk);
%     else
%         spk = [];
%         return;
%     end
% else
%     if(isfield(spk,'spiketimes') && ~isfield(spk,'spiketrain'))
%         spk = wz_get_SPKobj(spk.spiketimes);
%     end
end

%% define the spike kernel
krnl = wz_spk_kernel(ktype, kwdth);

% if spike times have been clipped, set spiketrain to NaN after last spike
% occurrence, but correct for the kernel width.
if(doclip == 1 && ~isempty(spk.clip))
    [spk.spiketimes, spk.trial_order] = wz_spk_clip(spk.spiketimes, spk.clip);

    for(t=1:spk.nTrials)
        %p = ~isnan(spk.spiketimes(t,:));
        spk.spiketimes(t, 1:sum(~isnan(spk.spiketimes(t,:)))) = spk.spiketimes(t, ~isnan(spk.spiketimes(t,:)));
%         spk.spiketrain(t,:)  =  histc(spk.spiketimes(t,~isnan(spk.spiketimes(t,:))), spk.xtime);
    end
    spk.clipped = true;
end

if(~isempty(spk))
    spkt = spk.spiketrain;
else
    spkt = [];
end

%% apply convolution on the spike train matrix
if(~isempty(spk) && ~isempty(spk.timewindow))
    spk.xtime = spk.timewindow(1) : spk.timewindow(2)+1;
    
    %spk.spikedensities = convn(spkt, krnl, 'same') .* 1000; % multiply with 1000 to get it from spikes/ms to spikes/s.
    
    for(i=1:size(spkt,1))
        spk.spikedensities(i,:) = conv(spkt(i,:), krnl, 'same') .* 1000; % multiply with 1000 to get it from spikes/ms to spikes/s.
    end    
    
    if(spk.clipped)
        spk.spikedensities(bsxfun(@gt, spk.xtime, spk.clip(:)+kwdth)) = NaN;
    end
    
    spk.Ndense = sum(~isnan(spk.spikedensities));
else
    spk.xtime = [];
    spk.spikedensities = [];
    spk.Ndense = [];
end

spk.kerneltype  = ktype;
spk.kernelwidth = kwdth;
spk.kernel      = krnl;
spk.nBoot = rep;
num_trial = size(spkt,1);  % might need correction to account for trials with different duration

%% get mean density
if(num_trial > 1)
    if(rep == 1)
        spk.meandensity   = nanmean(spk.spikedensities,1);
        spk.stedensity    = nanstd(spk.spikedensities,0,1)/sqrt(num_trial);

        spk.mediandensity = nanmedian(spk.spikedensities,1);
        spk.upperdense    = prctile(spk.spikedensities,75,1);
        spk.lowerdense    = prctile(spk.spikedensities,25,1);
        spk.IQRerr = 1.58 * (spk.upperdense - spk.lowerdense)  ./ sqrt(num_trial);
    else
        % use a bootstrap approach to determine mean and confidential intervals
        rng(sum(100*clock),'twister');

        bootmean     = nan(rep, size(spk.spikedensities,2));
        ctrials      = randi(num_trial,num_trial,rep); % wz_sample(trialvec, num_trial, 1, rep);
        ctrials(:,1) = 1:num_trial;

        for(b = 1:rep)
            cset = spk.spikedensities(ctrials(:,b),:);
            bootmean(b,:) = nanmean(cset,1);
        end

        spk.meandensity   = nanmean(bootmean,1);
        spk.stedensity    = nanstd(bootmean,0,1);

        spk.mediandensity = nanmedian(bootmean,1);
        spk.upperdense    = prctile(bootmean,75,1);
        spk.lowerdense    = prctile(bootmean,25,1);
        spk.IQRerr        = spk.upperdense - spk.lowerdense;
    end
else
    if(isempty(spk.spikedensities))
        spk.meandensity   = [];
        spk.stedensity    = [];
        spk.mediandensity = [];
        spk.upperdense    = [];
        spk.lowerdense    = [];
        spk.IQRerr        = [];
    else
        spk.meandensity   = spk.spikedensities;
        spk.stedensity    = zeros(size(spk.spikedensities));
        spk.mediandensity = spk.spikedensities;
        spk.upperdense    = zeros(size(spk.spikedensities));
        spk.lowerdense    = zeros(size(spk.spikedensities));
        spk.IQRerr        = zeros(size(spk.spikedensities));
    end
end


