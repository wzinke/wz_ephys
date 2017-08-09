function [coh, f] = pair_coh(csdmat, i, j, Fs);
%  
% ordinary coherence
% Squared Pairwise Coherence(the ith and jth channel) calculation 
% using ensemble average based on Multitaper Method
%
% Usage: 
%   [coh, f] = pair_coh(csdmat, i, j, Fs);
%    
% Input:
%   csdmat: spectrum density matrix, chans x chans x half_pad
%   i:   the ith channel selected
%   j:   the jth channel selected
%   Fs:  sampling frequency
% Output:
%    coh:  squared multiple coherence between ith channel and others(freq x 1)
%      f:  freq
% See also:
%   minor, 
%

%
%  Hualou Liang, 12/01/98, FAU
%  revised 7/14/99, change complex coh to real coh to save space for boot


[chans, chans, half_pad] = size(csdmat);


if nargin<4,
  Fs = 200;    % deault sampling rate
end

cohtmp = abs(csdmat(i, j, :)).^2 ./ (csdmat(i, i, :).*csdmat(j, j, :)); 
% coh = squeeze(cohtmp);
%% to make it real number to save space, 7/14/99
coh = abs(squeeze(cohtmp));

f = (0:half_pad-1)*Fs/(2*half_pad);  

% coh = abs(Pxy).^2 ./ (Pxx.*Pyy);
    															



