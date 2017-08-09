function  [kern,hcoeffs,hcoeffs2D] = hcoeff(peristim,rawdata,mfr,nstims,rsp,fast)


% *********************************************************************************************
% *********************************************************************************************
% *******   This code was authored by Michael Hill (buffalohill@me.com), it may not     *******
% *******   be used as a whole or in parts without the written permission of Michael    *******
% *******   Hill and any use of this code (in whole or in parts) must reference the     *******
% *******   paper in which the h-coefficient method was first introduced (contact me    *******
% *******   to inquire about appropriate referencing formats for this code).             ******
% *********************************************************************************************
% *********************************************************************************************

% bsl = [-1000,-200];
% rsp = [200,1000];
% nstims      =   1;
% mfr         =   3;
% fast        =   1;      	% Do the bootstrapped kernel method fast
%
% frombsl = round((bsl(1,1)/10) +500 + 1);
% tobsl = round((bsl(1,2)/10) +500);


dbstop if error


%% Generate the shuffled responses

rspdur = rsp(1,2)-rsp(1,1);                                                 % Response duration in ms
bslshuff = randshuff(nstims,rspdur,rawdata,mfr,fast);                       % Generate your shuffled responses

%% Generate a smooth PSTH from your response

corrkern = ones(nstims,1)*NaN;

for j = 1:nstims                                                            % Prep for the adaptive kernel analysis
    if j == 1
        forkern = peristim{1,1};
    else
        forkern = cat(2, forkern, peristim{j,1});
    end
    corrkern(j,1) = length(peristim{j,1});                                  % Preparation of the factor by which you need to multipy the ssvkern output, as it is normalized to the mean firing rate
end

[kern(2,:),kern(1,:)] = ssvkernel(forkern.*10^(-3),(-4.99:0.01:5));         % Smooth your PSTH

kern(2,:) = kern(2,:) * (squeeze(nanmean(corrkern)));                       % Reverse the normaliztion applied by ssvkernel and transform to Hz.

if isempty(kern(2,~isnan(kern(2,:))))                                       % Check for the presence of any firing rate whatsoever
    kern(2,:) = zeros(1,length(kern));
end

%% Find your local extrema

fromrsp = round((rsp(1,1)/10) +500 + 1);                                    % Find the index of the local max and local min of the response for further analysis
torsp = round((rsp(1,2)/10) +500);

indmax = find(kern(2,fromrsp:torsp) == max(kern(2,fromrsp:torsp)),1,'first') + (fromrsp - 1);
indmin = find(kern(2,fromrsp:torsp) == min(kern(2,fromrsp:torsp)),1,'first') + (fromrsp - 1);

if isempty(indmax) || isempty(indmin)
    indmax = NaN;
    indmin = NaN;
end

%% Calculate the h-coefficient

[hcoeff_max,hcoeff_min,hcoeff_max2D,hcoeff_min2D] = hcoeff_gen(bslshuff,kern(2,:),indmin,indmax,mfr);
hcoeffs = [hcoeff_max,hcoeff_min];
hcoeffs2D = [hcoeff_max2D;hcoeff_min2D];

end
function [bslshuff] = randshuff(nstims,rspdur,rawdata,mfr,fast,nshuff)
%
% randshuff is a function that takes all spikes from within one whole
% recording session (rawdata) and from it creates nshuff (i.e. usually 1000
% due to time constraints) randomly generated responses irrespective of whether
% there actually was a real response during some part or not. To each one
% of these shuffled responses randshuff applies adptive kernel
% smoothing (ssvkernel) and then calculates a master response M across all
% of them, from which one can then later calculate the h-coefficient.
% To execute ssvkernel 1000 times takes a very long time on a normal PC, so
% randshuff also offers a 'fast' method (fast = 1). Instead of cutting out
% 1000 segments of data of a certain length and then applying ssvkernel to
% each segment, the fast method cuts out 100 segments of 10 times that
% length. It then only needs to apply ssvkernel 100 times instead of 1000
% times. After running ssvkernel 10 times randshuff then cuts each smoothed
% segment into 10 parts, which again results in 1000 segments. Running
% randshuff in fast mode is not as good as running it in 'slow' (fast = 0)
% mode but in most cases it's ok as long as your interstimulus intervalls
% are at least 10 times longer then your baseline periods.
%
% Input:
%
% nstims:
% nstims is the number of trials or stimuli that were run during the whole
% experiments
%
% rspdur:
% rspdur is the response duration in ms, i.e. the duration of the individual
% segments that will make up your bootstraped responses
%
% rawdata
% rawdata is a vector containng timestamps for each spike that was recorded
% during yoru experiment.
%
% mfr
% mfr is the mean firing rate for this neuron across the whole
% experiment
%
% fast
% fast = 1 runs ranshuff in fast mode. fast = 0 runs randshuff in slow mode

dbstop if error
if nargin < 5
    fast = 1;
end
if nargin < 6
    nshuff = 1000;
end
kernshuff = randperm(nshuff);
shuffsy = zeros(nstims,1) * NaN;
nsegs = floor(((floor(rawdata(end))-2000)-(ceil(rawdata(1,1))))/(rspdur*10));      % Add 2 x 1000 = 2000 ms here so that you can always subtract them again afterwards to get rid of any edge artifact in the ssvkernel
if fast
    for i = 1:nshuff/10
        clc
        disp([num2str(((i-1)*10)+1),' to ',num2str(i*10),' out of ',num2str(nshuff),' shuffles being processed with the "fast" method...'])
        permind1 = randperm(nsegs*10);
        permind2 = randperm(nsegs);
        for j = 1:nstims
            shuffsy(j,1) = length(rawdata((rawdata > ((permind1(1,j)-1)*rspdur) & (rawdata <= (((permind1(1,j)-1)*rspdur)+rspdur)))));
            if j == 1 % Prep for the adaptive kernel analysis
                forkern = (rawdata((rawdata > ((permind2(1,j)-1)*(rspdur*10))) & (rawdata <= ((permind2(1,j)-1)*(rspdur*10))+(rspdur*10)+2000)) - ((permind2(1,j)-1)*(rspdur*10)));
            else  % Only do this 100 times, but with a x10 longer duration to save time (instead of 100 times)
                forkern = cat(2,forkern,(rawdata((rawdata > ((permind2(1,j)-1)*(rspdur*10))) & (rawdata <= ((permind2(1,j)-1)*(rspdur*10))+(rspdur*10)+2000)) - ((permind2(1,j)-1)*(rspdur*10))));
            end
        end
        [prekern{i}(2,:),prekern{i}(1,:)] = ssvkernel(forkern.*10^(-3),(0.01:0.01:(((rspdur*10)+2000)/1000))); % Execute your adaptive kernel analysis
        prekern{i}(2,:) = prekern{i}(2,:) .* 10 .* (squeeze(nanmean(shuffsy))*(1000/rspdur));
        for j = 1:10
            k = ((i-1)*10)+j;
            bslshuff.kern{kernshuff(k)}(1,:) = (0.01:0.01:rspdur/1000);
            bslshuff.kern{kernshuff(k)}(2,:) = prekern{i}(2,(j-1)*(rspdur/10)+1:j*(rspdur/10));      % Reverse the normaliztion applied by ssvkernel and transform to Hz. Also cut off the first and last 50 ms, as these have edge artifact. You can do so because above you've added 100 ms to the baseline, so youre going to end up with the right length again.
        end
    end
else
    for i = 1:nshuff
        clc
        disp([num2str(i),' out of ',num2str(nshuff),' shuffles being processed with the "slow" method (this might take a while)...'])
        permind = randperm(nsegs*10);
        for j = 1:nstims
            shuffsy(j,1) = length(rawdata((rawdata > ((permind(1,j)-1)*rspdur) & (rawdata <= (((permind(1,j)-1)*rspdur)+rspdur)))));
            if j == 1   % Prep for the adaptive kernel analysis
                forkern = (rawdata((rawdata > ((permind(1,j)-1)*rspdur)) & (rawdata <= ((permind(1,j)-1)*rspdur)+rspdur+2000)) - ((permind(1,j)-1)*rspdur))';
            else
                forkern = cat(2,forkern,(rawdata((rawdata > ((permind(1,j)-1)*rspdur)) & (rawdata <= ((permind(1,j)-1)*rspdur)+rspdur+2000)) - ((permind(1,j)-1)*rspdur))');
            end
        end
        [prekern{i}(2,:),prekern{i}(1,:)] = ssvkernel(forkern.*10^(-3),(0.01:0.01:((rspdur+2000)/1000))); % Execute your adaptive kernel analysis
        prekern{i}(2,:) = prekern{i}(2,:) .* 10 .* (squeeze(nanmean(shuffsy))*(1000/rspdur));
        bslshuff.kern{kernshuff(i)}(1,:) = (0.01:0.01:rspdur/1000);
        bslshuff.kern{kernshuff(i)}(2,:) = prekern{i}(2,:);
    end
end
clc

%% Now go through your shuffled adaptively smoothed baselines/responses and look for responses
bslshuffmaxes = zeros(10001,nshuff);
bslshuffminses = bslshuffmaxes;
for a = 1:nshuff  % for each of your shuffled samples
    kern = (bslshuff.kern{a}(2,:)./mfr)-1;                                  % Normalize it to your channel's mean firing rate and subtract 1 so that your baseline value is = 0
    indmax = find(kern == max(kern),1,'first');                             % Find the areas for your local maximum and save them into the maxes variable.
    if ~isempty(indmax) && kern(1,indmax) > 0                               % just to make sure your maximum is above the mean firing rate of this channel, as else this analysis makes no sense
        firstdec = norm(floor(norm(kern(1,indmax))*10)/10);                 % Find out what the next decimal point is between your max and the channel baseline (in units of channel baseline firing rate, e.g. if the max firing rate is 1.568 times the mean firing rate of the channel, then your next decimal point towards that mean firing rate is 1.5)
        u = indmax(1); v = u;
        for b = 1:round(firstdec*10)+1
            while u > 1 && kern(1,u) > roundn(firstdec-((b-1)*0.1),-1)      % This travels along the curve towards the left (heading downwards from the local maximum) until it hits the next decimal multiple of the mean firing rate of the whole channel
                u = u-1;
            end
            while v < length(kern) && kern(1,v) > roundn(firstdec-((b-1)*0.1),-1)   % This travels along the curve towards the right (heading downwards from the local maximum) until it hits the next decimal multiple of the mean firing rate of the whole channel
                v = v+1;
            end
            bslshuffmaxes(end-(round(firstdec*10)+1)+b,a) = trapz((kern(1,u:v))-roundn(firstdec-((b-1)*0.1),-1)).*0.01; % "-firstdec" makes sure the only the area between the next decimal multiple of the baseline firing rate and the peack that reaches just above this threshold is calculated. The factor 0.1 is there because the spacing of our x-axis is 0.1, see help trapz for more on this
        end
    end
    indmin = find(kern == min(kern),1,'first');                             % Find the areas for your local minimum and save them into the minses variable.
    if ~isempty(indmax) && kern(1,indmin) < 0                               % just to make sure your maximum is above the mean firing rate of this channel, as else this analysis makes no sense
        firstdec = norm(floor(norm(kern(1,indmin))*10)/10);                 % Find out what the next decimal point is between your max and the channel baseline (in units of channel baseline firing rate, e.g. if the max firing rate is 1.568 times the mean firing rate of the channel, then your next decimal point towards that mean firing rate is 1.5)
        u = indmin(1); v = u;
        for b = 1:round(firstdec*10)+1
            while u > 1 && kern(1,u) < -roundn(firstdec-((b-1)*0.1),-1)     % This travels along the curve towards the left (heading upwards from the local minimum) until it hits the next decimal multiple of the mean firing rate of the whole channel
                u = u-1;
            end
            while v < length(kern) && kern(1,v) < -roundn(firstdec-((b-1)*0.1),-1) % This travels along the curve towards the right (heading upwards from the local minimum) until it hits the next decimal multiple of the mean firing rate of the whole channel
                v = v+1;
            end
            bslshuffminses(end-(round(firstdec*10)+1)+b,a) = norm(trapz(kern(1,u:v) + roundn(firstdec-((b-1)*0.1),-1)).*0.01); % "+(firstdec-((b-1)*0.1))" makes sure that only the area between the next decimal multiple of the baseline firing rate and the peack that reaches just below this threshold is calculated. The factor 0.1 is there because the spacing of our x-axis is 0.1, see help trapz for more on this
        end
    end
end
bslshuffmaxes = bslshuffmaxes(2:end,:) - bslshuffmaxes(1:end-1,:);             % This subtracts the row above from each row. Before doing so your area measurements represent the area above each decima multiple of the baseline firing rate. After doing this, your area measurements are now corrected to reflect only the area inbetween this and the next following decimal. In other words, instead of reporting how large the area is above e.g. 1.7 * meanfr(channel), we now are reporting how large the area is inbetween 1.7 * meanfr(channel) and 1.8 * meanfr(channel), i.e. the maximum in each band
bslshuffminses = bslshuffminses(2:end,:) - bslshuffminses(1:end-1,:);
bslshuff.maxes = max(bslshuffmaxes,[],2);
bslshuff.minses = max(bslshuffminses,[],2);
end
function [hcoefficient_max,hcoefficient_min,hcoefficient_max2D,hcoefficient_min2D] = hcoeff_gen(bslshuff,kern,indmin,indmax,normit)
maxes = zeros(10001,1);                                                     % maxes & minses are two long vector of zeros, into which we can then enter the area measurements during the h-coefficient calculation further down
minses = maxes;
kern = (kern/normit)-1;                                                     % Normalize the response to the channel's baseline value and subtracts 1 so that the baseline value is = 0
if ~isnan(indmax) && kern(1,indmax) > 0                                     % This makes sure the maximum is above the mean firing rate of this channel, as else this analysis makes no sense
    firstdec = norm(floor(norm(kern(1,indmax))*10)/10);                     % Finds out what the next decimal point is between the local max and the channel baseline (in units of channel baseline firing rate, e.g. if the max firing rate is 1.568 times the mean firing rate of the channel, then your next decimal point towards that mean firing rate is 1.5)
    u = indmax(1); v = u;
    for b = 1:round(firstdec*10)+1
        while u > 1 && kern(1,u) > roundn(firstdec-((b-1)*0.1),-1)          % This travels along the curve towards the left (heading downwards from the local maximum) until it hits the next decimal multiple of the mean firing rate of the whole channel
            u = u-1;
        end
        while v < length(kern) && kern(1,v) > roundn(firstdec-((b-1)*0.1),-1) % This travels along the curve towards the right (heading downwards from the local maximum) until it hits the next decimal multiple of the mean firing rate of the whole channel
            v = v+1;
        end
        maxes(end-(round(firstdec*10)+1)+b,1) = trapz((kern(1,u:v))-roundn(firstdec-((b-1)*0.1),-1)).*0.1; % "-firstdec" makes sure that only the area between the next decimal multiple of the baseline firing rate and the peack that reaches just above this threshold is calculated. The factor 0.1 is there because the spacing of our x-axis is 0.1, see help trapz for more on this
    end
end
if ~isnan(indmin) && kern(1,indmin) < 0                                     % This makes sure the minimum is below the mean firing rate of this channel, as else this analysis makes no sense
    firstdec = norm(floor(norm(kern(1,indmin))*10)/10);                     % Finds out what the next decimal point is between your min and the channel baseline (in units of channel baseline firing rate, e.g. if the max firing rate is -0.568 times the mean firing rate of the channel, then your next decimal point towards that mean firing rate is 0.5)
    u = indmin(1); v = u;
    for b = 1:round(firstdec*10)+1
        while u > 1 && kern(1,u) < -roundn(firstdec-((b-1)*0.1),-1)         % This travels along the curve towards the left (heading upwards from the local minimum) until it hits the next decimal multiple of the mean firing rate of the whole channel
            u = u-1;
        end
        while v < length(kern) && kern(1,v) < -roundn(firstdec-((b-1)*0.1),-1) % This travels along the curve towards the right (heading upwards from the local minimum) until it hits the next decimal multiple of the mean firing rate of the whole channel
            v = v+1;
        end
        minses(end-(round(firstdec*10)+1)+b,1) = norm(trapz(kern(1,u:v)+roundn(firstdec-((b-1)*0.1),-1)).*0.1); % "+(firstdec-((b-1)*0.1))" makes sure that only the area between the next decimal multiple of the baseline firing rate and the peack that reaches just below this threshold is calculated. The factor 0.1 is there because the spacing of our x-axis is 0.1, see help trapz for more on this
    end
end
maxes = maxes(2:end,1) - maxes(1:end-1,:);                                  % From each row this subtracts the row right above it. Before doing so the area measurements represent the area above each decimal multiple of the baseline firing rate. After doing this, the area measurements are now corrected to reflect only the area inbetween this and the next following decimal. In other words, instead of reporting how large the area is above e.g. 1.7 * meanfr(channel), we now are reporting how large the area is inbetween 1.7 * meanfr(channel) and 1.8 * meanfr(channel), i.e. in each band
minses = minses(2:end,1) - minses(1:end-1,:);
hcoefficient_max = sum(maxes > bslshuff.maxes)/sum(bslshuff.maxes > 0);     % This now first counts in how many bands the response shows a bigger area (i.e. is wider) then the max in bands across the shuffstrapped samples. It then divides this number by the number of bands that showed any response at all in the shuffstrapped samples (i.e. the maximum reach of the highest baseline response). The resulting number means the following: If it is smaller then 1, then this means that the response was wider then the baseline responses in only some of the bands (although it may overall reach higher then any other, in which case hcoefficient < 1 would mean that at it's base this high response waas still narrower then the maxes from the baselines). If hcoefficient = 1 this means that the response was the same hight as the highest response in the baseline but wider in every band. If hcoefficient > 1 this means that the response was both higher then the baseline reseponses and wider in more bands.Any hcoefficient above zero can be reagegarded as significant dependent on your criteria. hcoefficient is higher then one should definitively be regarded as highly significant
hcoefficient_min = sum(minses > bslshuff.minses)/sum(bslshuff.minses > 0);
prov = (nnz(maxes)-nnz(bslshuff.maxes));                                    % number of bands that reac higher then in the shuffstrapped baseline
if prov > 0
    hcoefficient_max2D(1,1) = prov;
else
    hcoefficient_max2D(1,1) = 0;
end
prov = (nnz(maxes)-nnz(bslshuff.maxes));                                    % number of bands that reac higher then in the shuffstrapped baseline
if prov > 0
    hcoefficient_max2D(1,1) = prov;
else
    hcoefficient_max2D(1,1) = 0;
end
hcoefficient_max2D(1,2) = sum(maxes > bslshuff.maxes) - hcoefficient_max2D(1,1);  % number of bands that aren't higher but are wider then the in the shuffled baseline
hcoefficient_max2D(1,3) = sum(bslshuff.maxes > 0);                          % number of bands in the shuffled baseline
prov = (nnz(minses)-nnz(bslshuff.minses));                                  % number of bands that reach higher then in the shuffled baseline
if prov > 0
    hcoefficient_min2D(1,1) = prov;
else
    hcoefficient_min2D(1,1) = 0;
end
hcoefficient_min2D(1,2) = sum(minses > bslshuff.minses) - hcoefficient_min2D(1,1);  % number of bands that aren't higher but are wider then the in the shuffled baseline
hcoefficient_min2D(1,3) = sum(bslshuff.minses > 0);                         % number of bands in the shuffled baseline
end

