function [coh, f] = part_coh(csdmat, i, j, Fs);
%  
% Squared Partial Coherence(the ith and jth channel) calculation 
% using ensemble average based on Multitaper Method
% For details, see Jenkins and Watts(1968), p. 489
%
% Usage: 
%   [coh, f] = part_coh(csdmat, i, j, Fs);
%    
% Input:
%   csdmat: spectrum density matrix, chans x chans x half_pad
%   i:   the ith channel selected
%   j:   the jth channel selected
%   Fs:  sampling frequency
% Output:
%    coh:  squared multiple coherence between ith channel and others
%          (1 x freq) complex vector
%      f:  freq
% See also:
%   minor, 
%

%
%  Hualou Liang, 11/09/98, FAU
%  Revised 01/06/99, remove waitbar for large data set
%


[chans, chans, half_pad] = size(csdmat);


if nargin<4,
  Fs = 200;    % deault sampling rate
end


coh = [];
Mi =  minor(csdmat, i, i);
Mj =  minor(csdmat, j, j);
Mij =  minor(csdmat, i, j);
% h = waitbar(0,'Please wait...');
for k=1:half_pad,
  % tmp = det(Mij(:,:,k)) / sqrt(det(Mi(:,:,k)) * det(Mj(:,:,k)));
  tmp = det(Mij(:,:,k))^2 / (det(Mi(:,:,k)) * det(Mj(:,:,k)));
  coh = [coh tmp];
%  waitbar(k/half_pad)
end
% close(h);

f = (0:half_pad-1)*Fs/(2*half_pad);  

% coh = abs(Pxy).^2 ./ (Pxx.*Pyy);
    															




