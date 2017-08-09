function power = autopower(aredat,arndat,n,fs)
% AUTOPOWER Compute the auto power
% 
% Usage: 
%   power = autopower(aredat,arndat,n,fs);
% 
%   aredat: AR coefficient file
%   arndat: AR noise file
%   n: Number of frequency bins
%   fs: Sampling rate
% 
% Output(s):
%   power - auto power in [Time_points x Frequency_bins x Channels]

% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.3$ $Date: 13-Sep-2007 17:04:26$
% SHIS UT-Houston, Houston, TX 77030, USA.
%
% Hualou Liang

csd4d = MAR_csd4dmem(aredat,arndat,n,fs);
[paircoh, partcoh, mulcoh, autospect] = MAR_coh3D(csd4d, fs);
% change into power
power = abs(autospect).^2;

end%autopower

% [EOF]