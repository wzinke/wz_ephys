function [mar] = MAR_make(coeff,noise);
%
% To make it more general for any channel numbers
%
% A batch version, to covert MAR output to MAR structure
%
% Usage:
%    [mar] = MAR_make(coeff,noise);
% Example: 
%    [mar] = MAR_make(coeff,noise);
%    [x] = ensembding(mar, T, chans, N)
% Input: 
%    coeff: row vector(1 x 1350)
%    noise: row vector (1 x 225)
% MAR is modeled as 
%   A0*X_t + A1*X_(t-1) + ... + A5*X_(t-5) = Et 
%
% Output:
%    mar: structure
% 
% See also MAR_INIT_15
%

% 
% Hualou Liang, 01/26/99, FAU
% revised 3/19/99
%


% load dat
order = length(coeff)/length(noise) - 1;
chans = sqrt(length(noise));

A = [];
for i=1:order,
  % for i=1:5,
  % A = [A reshape(coeff(i*225+1:(i+1)*225),15,15)'];
  A = [A reshape(coeff(i*chans*chans+1:(i+1)*chans*chans),chans,chans)']; 
end;

% Reshaped data format ONLY for arsim.m, which mean A = - A returned
%  in order to fit the format requirement of arsim
% A = -A;   % for prepare for arsim

% Noise coefficient matrix
% C = reshape(noise,15,15)';
C = reshape(noise,chans,chans)';

% generate structure MAR
% mar = mar_init_15(A, C);
mar = mar_init(A, C);

