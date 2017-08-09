
function [paircoh, partcoh, mulcoh, autospect] = ...
	MAR_coh3D(csd4d, Fs);

% The only difference with coh3D is without global variables and t, f 
% calculate EVOLUTIONARY ORDINARY, PARTIAL, MULTIPLE COHERENCE & AUTOSPECTRUM 
% based on csd4D spectral matrix
% 
% Usage:
%  [paircoh, partcoh, mulcoh, autospect] =MAR_coh3D(csd4d, Fs); 
% Example:
%   make sure to set global variables before running it, see tfcohmain.m
%   [paircoh, partcoh, mulcoh, spect, t, f] = coh3D(csd4d,t,f);
% Input:
%   csd4d: spectral matrix, t x chans x chans x f
%   Frng:    freq range(Hz)
%   Trng:    time range(msec)
%   N: the numbers of trial(in case 887)
% Output:
%   paircoh:  ordinary coherence matrix, t x f x CombinedChans
%   partcoh:  partial coherence
%   mulcoh:   multiple coherence
%   autospect:  auto spectral of each channel
%   Frng:    freq range(Hz)
%   Trng:    time range(msec)
%

%  
%  Hualou Liang, 01/07/99, FAU
% 




[TT, chans, chans, FF] = size(csd4d);
combchan = chans*(chans-1)/2; 

paircoh = zeros(TT, FF, combchan); % t x f x combchan
partcoh = zeros(TT, FF, combchan); % t x f x combchan
mulcoh = zeros(TT, FF, chans); % t x f x chans
autospect = zeros(TT, FF, chans); % t x f x chans

%h = waitbar(0,'Please wait...');
for tt=1:TT,
  tt;
  csd = squeeze(csd4d(tt,:,:,:));
  
  cohtmp1 = []; 
  cohtmp2 = []; 
  cohtmp3 = []; 
  specttmp = [];  % for autospectra

  %%% calculate ordinary and partial coherence

  for i=1:chans,
	for j=i+1:chans,
	  [coh1, f] = pair_coh(csd,i,j,Fs);  % TODO: need new func to remove f
	  [coh2, f] = part_coh(csd,i,j,Fs);  % TODO: need new func to remove f
	  cohtmp1 = [cohtmp1 coh1];
	  cohtmp2 = [cohtmp2 coh2'];
	end
  end
  paircoh(tt,:,:) = cohtmp1;
  partcoh(tt,:,:) = cohtmp2;

 %%% calculate multiple coherence & autospectrum 
  for i=1:chans,
	[coh3, f] = mul_coh(csd, i, Fs);
	cohtmp3 = [cohtmp3 coh3'];
  end
  mulcoh(tt,:,:) = cohtmp3;

 %%% calculate autospectrum 
  for i=1:chans,
	specttmp = [specttmp squeeze(csd(i,i,:))];
  end
  autospect(tt,:,:) = specttmp;

  %waitbar(tt/TT)
end
%close(h);
  
