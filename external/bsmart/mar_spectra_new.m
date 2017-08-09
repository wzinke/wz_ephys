function [pow,f, H] = mar_spectra_new (mar, N, Fs);

% The only difference between current version and mar_spectra.m is
% the output of dim chan x chan x freq
%  old version: freq x chan x chan
% Transfer function also calculated 
%
% This change is for calculation of multiple/parital coherence and the like
%
%  [pow,f, H] = mar_spectra_new (mar, N, Fs);
%  
% Get Hermitian PSD from MAR coefficients
% See Marple (1987; page 408)
%  Xt + A1*X_(t-1) + ... + = C*W_t
%
% Input:
%   mar.lag(k).a:     is ar coefficient matrix at lag k
%   mar.noise_cov:    estimated noise deviation
%   Fs:    samples per second
%   N: length of returned spectrum
% Output:
%   pow:  spectra matrix, dim x dim x freq (here is different from old one)
%   H:   transfer function(square), dim x dim x freq
%   f: freqency
%

%
% Hualou Liang, 01/26/99, FAU
%


order = mar.order;
d=size(mar.lag(1).a,1);
chans = size(mar.noise_cov,1);

ff=0;
for w=pi/N:pi/N:pi,
  ff=ff+1;
  af_tmp=eye(d);
  for k=1:order,
    af_tmp=af_tmp+mar.lag(k).a*exp(-i*w*k);
  end
  iaf_tmp=inv(af_tmp);  % 3 x 3 
  h_tmp = abs(iaf_tmp).^2;
  H(:,:,ff) = h_tmp ./ (sum(h_tmp,2)*ones(1,chans)); % transfer function
  pow(:,:,ff) = iaf_tmp * mar.noise_cov * iaf_tmp';
end

% f=[0.5*Fs/N:0.5*Fs/N:0.5*Fs]'; 
f=[0:N-1]*Fs/(2*N);
  

