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

if(exist('limit','var') == 0 || isempty(limit) == 1)
    limit = 150;
end

if(exist('timwin','var') && ~isempty(timwin))
    spktimes(spktimes>timwin(2) | spktimes < timwin(1)) = NaN;
end
num_trial = size(spktimes,1);

corr_times = [];
bins = [(limit*(-1)) : limit];

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

