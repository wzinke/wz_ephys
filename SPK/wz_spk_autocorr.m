function [out_hist, xval, corr_times] = wz_spk_RACH(spktimes, timwin, limit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on the RCCH function written by Alex Thiele
% 
%    calculate auto-correlation histogram based on trial-wise spike times.
%    this utilizes the “raw auto-coincidence histogram” (RACH).
%
%
% modified by wolf zinke


%% check data input
if(exist('timwin','var') == 0)
    timwin = [0;600];
elseif(isempty(timwin) == 1)
    timwin = [0;600];
end

if(exist('limit','var') == 0)
    limit = 150;
elseif(isempty(limit) == 1)
    limit = 150;
end

num_trial = size(spktimes,1);

spktimes(spktimes>timwin(2) | spktimes < timwin(1)) = NaN;
% spkcnt = sum(isfinite(sptspk),2);

corr_times = [];
bins = [(limit*(-1)) : limit];

spktrain = zeros(num_trial,diff(timwin)+1);
spk_x    = timwin(1) : timwin(2);

actrial = nan(num_trial,2*limit+1);

for (i=1: num_trial)
    x = spktimes(i, ~isnan(spktimes(i,:)));
    
    xh = histc(x, spk_x);
    
    actrial(i,:) = xcorr(xh, limit);
   
end

%ac = fliplr(conv(spktrain(i,:),fliplr(spktrain(i,:))))

for (i=1: num_trial)
   x = spktimes(i, ~isnan(spktimes(i,:)));

   if (~isempty(x))
      X = repmat(x',1,size(x,2));
      Y = repmat(x,   size(x,2),1);
      Z = Y - X;
      Z = Z(:);
      corr_times = [corr_times Z( (Z >= limit * (-1)) & (Z <= limit) )'];
   end
end

corr_times(corr_times == 0) = [];  % Get rid of these expected values
[out_hist, xval] = hist(corr_times,bins);

disp('done');
