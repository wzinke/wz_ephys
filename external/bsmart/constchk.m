function [R1, R2, diffR, ratio] = constchk(dat1, dat2);
% function [R, Ru, Rd, trueR, denoRm] = consistncycheck1(fname, mu);

%%% Using all one pool of size 888
%
% For each data window(16), computing all chan-pairwise cross
% correlation function at different time lag(NTOT), then check each value of 
% cross-correlation whether it's higer than certain signicant threshold(NSIG), 
% the ratio of NSIG to NTOT is probability 
% Usage:
%   
% Input:
%   fname1: onewin.dat (3D data), monkey dat
%   fname2: smimulation data based on MAR, same format as fname1
%   mu: 1 x len (wrong: 11 x 120)
% Output: ( are scalar number)
%   R1: radius of monkey dat
%   R2: radius of simulation dat
%   diffR: difference between R1 and R2
%   ratio:  diffR/R1
%

% 
% Hualou Liang, 03/16/99, FAU
%



% [dat1] = readchkdat(fname1);
% [dat2] = readchkdat(fname2);

[val1] = corrlag(dat1);  % 1 x 1320
[val2] = corrlag(dat2);

  

% estimate distance from each row to origin(0)
R1 = sqrt(sum(val1.^2));  % 1 x 1, for monkey
R2 = sqrt(sum(val2.^2));  % 1 x 1

diffR = sqrt(gtm_dist(val1, val2));  % 1 x 1

ratio = diffR/R1;

% trueR = sqrt(sum(mu.^2));  % a scalar number


% stdX = std(prob)/sqrt(N);  % stderror


%function [dat] = readchkdat(fname);

%ntrls = 314;
%NPTS = 16;
%NCHN = 15;

%fid=fopen(fname,'rb','ieee-le');
%dat=fread(fid, ntrls*NCHN*NPTS, 'float');
%dat = reshape(dat,[NPTS NCHN ntrls]);
%fclose(fid);


function [val] = corrlag(dat);
global WIN NPTS NTRLS NCHN 
% calculate correlation for various lags, and put in one row vector
% Input:
%  dat: 16 x 15 x 888
% Output:
%  val: 1 x Len
%


% runs = 1;
%NCHN = 15;
B=NTRLS;     % size of pool
nlags = 5;

sel = nonzeros(tril(reshape(1:NCHN*NCHN, NCHN, NCHN))); % select

% disp('Bootstrap 100 times, pool size of 200');
  % idx = randperm(ntrls);
  % dat1 = dat(:,:,idx(1:B));
  val = 0;
  for j=1:B,
%	val = val + xcorr(dat1(:,:,j), nlags,'coeff');  % (2*nlags-1) x 15^2
	val = val + xcorr(dat(:,:,j), nlags,'coeff');  % (2*nlags-1) x 15^2
  end  
  val = val/B;
  val = val(:,sel);
  val = val(:)';     % runs (1) x p [No_of_variables(NoV)]  





