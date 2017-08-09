function [partsp] = part_spect(csd4d, i, j, varargin );

% Calculate partial auto/cross-spectr between channel i and j conditoned on 
% fixed channel k1, [k2]
% Usage:
%   [partsp] = part_spect(csd4d, i, j, k);
%   [partsp] = part_spect(csd4d, i, j, k1, k2,... );
% Example:
%   [partsp] = part_spect(csd4d, 2, 8, k1, k2);
%    2 - striate2, 5 - prestriate2, 8 - motor2 
% Input:
%   csd4d: spectral matrix, t x chans x chans x f
%   t:  time range in msec
%   f:  frequency rage in Hz
%   i,j,k:  channel number
% Output:
%   partsp: partial auto/cross-spectra between chan i and j 
%           conditioned on channel k, t x f 
%   t:  time range in msec
%   f:  frequency rage in Hz
% 

%
%  Hualou Liang, 02/07/99, FAU
% 


global Fs

% nout = max(nargout,1)-1;

% select channels i, j, k
if i==j,  % auto-spect
  csd4d = csd4d(:,[i varargin{:}], [i varargin{:}], :);
  [TT, chans, chans, FF] = size(csd4d);   % t x 3 x 3 x f
  partsp = zeros(TT, FF);  % t x f
  for t=1:TT,
	csd1 = squeeze(csd4d(t,:,:,:)); 
	csd = minor(csd1, 1, 1);
	for k=1:FF,
	  tmp = inv(csd(:,:,k));
	  partsp(t,k)=csd1(1,1,k) - csd1(1, 2:end, k)*tmp * csd1(2:end,1,k);
	end
  end
else      % cross-spec
  csd4d = csd4d(:,[i j varargin{:}], [i j varargin{:}], :);
  [TT, chans, chans, FF] = size(csd4d);   % t x 3 x 3 x f
  partsp = zeros(TT, FF);  % t x f
  for t=1:TT,
	csd = minor(squeeze(csd4d(t,:,:,:)), 2, 1);
	for k=1:FF,
	  tmp = inv(csd(2:end,2:end,k));
	  partsp(t,k)=csd(1, 1, k) - csd(1, 2:end,k)*tmp * csd(2:end,1,k);
	end
  end 
end








