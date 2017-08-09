function [coh, f] = mul_coh(csdmat, ith, Fs);
%  
% Squared Multiple Coherence(the ith chan and others) calculation 
% using ensemble average based on Multitaper Method
% For details, see Jenkins & Watts(1968), P.487
%
% Usage: 
%   [coh, f] = mul_coh(csdmat, ith, Fs);
%    
% Input:
%   csdmat: spectrum density matrix, chans x chans x half_pad
%   ith:   the ith channel selected
%     Fs:  sampling frequency
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


if nargin<3,
  Fs = 200;    % deault sampling rate
end

coh = [];
Mi =  minor(csdmat, ith, ith);
% h = waitbar(0,'Please wait...');
for i=1:half_pad,
  tmp = 1-det(csdmat(:,:,i)) / (csdmat(ith, ith, i) .* det(Mi(:,:,i)) );  
        % complex number
  coh = [coh tmp];
%  waitbar(i/half_pad)
end
% close(h);

f = (0:half_pad-1)*Fs/(2*half_pad);  

% coh = abs(Pxy).^2 ./ (Pxx.*Pyy);
    															





