function [csd4d, freq] = MAR_csd4dmem(arcoeff, arnoise, N, Fs);

% calculate spectral matrix, csd4d, based on output of MARfit 
%
% Usage:
%   [csd4d, H, freq] = MAR_csd4d(arcoeff, arnoise, N, Fs);
% Input:
%   arcoeff: t x Len, Len - length of one window MAR coefficients
%   arnoise: t x 225
%   N:  length of returned spectrum
%   Fs:  samples per second
%   time:   T x 1, time labels
% Output:
%   csd4d: spectral matrix, t x chans x chans x f
%   H: transfer function, t x chans x chans x f
%   time: time range(msec)
%   freq: freq range(Hz)
% 

%
% Hualou Liang, 01/26/99, FAU
% revised 01/28/99
%


[T, len] = size(arnoise);
chans = sqrt(len);

%csd4d = zeros(T, chans, chans, N);  % t x chans x chans x f
% H = zeros(T, chans, chans, N);  % t x chans x chans x f

% h = waitbar(0,'Please wait...');
for t=1:T,
  t;
  mar = MAR_make(arcoeff(t,:), arnoise(t,:));
  [csd, freq, Htmp] = mar_spectra_new (mar, N, Fs);
  csd4d(t,:,:,:) = csd;
%  H(t,:,:,:) = Htmp;
  
  % waitbar(t/T);
end
% close(h);



