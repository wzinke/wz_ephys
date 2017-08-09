function coherence = paircoherence(aredat,arndat,n,fs)
% PAIRCOHERENCE Compute the coherence pairs
% 
% Usage: 
%  [coherence] = paircoherence(aredat,arndat,n,fs);
% 
% Input(s):
%  aredat: AR coefficient file
%  arndat: AR noise file
%  n: Number of frequency bins
%  fs: Sampling rate
% 
% Output(s):
%   coherence   -   Pairwised ordinary coherence in the format of 
%                   [Time_points x Frequency_bins x Coherence]; the
%                   sequence of coherence follows: (1,1),(1,2),...,(1,k),
%                   (2,2),(2,3),...,(2,k),...,(k,k), where k is the number
%                   of channels, amounting to 1/2*k(k-1) coherence
%                   measures.
% 
% Example:
% 
% See also: .

% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.2$ $Date: 12-Sep-2007 11:35:38$
% SHIS UT-Houston, Houston, TX 77030, USA.
% 
% Lei Xu, Hualou Liang

csd4d     = MAR_csd4dmem(aredat,arndat,n,fs);
paircoh   = MAR_coh3D(csd4d,fs);
coherence = paircoh;

end%function

% [EOF]