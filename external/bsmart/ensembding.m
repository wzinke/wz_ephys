function [x] = ensembding(mar, T, chans, N)
%
% Generate ensemble dat of size N 
%
%  Usage:
%	dat = ensembding(mar, T, chans, N);
%  Input:
% 	N:  size of ensemble data set
% 	T:  length of dat x
%   chans: the number of channels
%  Output:
%    x: T x chans x N matrix
%
%  See also: 
%    mar_gen 
%

%
%  Hualou Liang, 09/10/98, FAU
%

x = zeros(T, chans, N);
%h = waitbar(0,'Please wait...');
for i=1:N,
   dat = mar_gen(mar, T);	% dat: 3xT matrix
   x(:,:, i) = dat';
   %waitbar(i/N)
end;
%close(h);

