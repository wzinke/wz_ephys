function [mar] = mar_init(A, C);

% To make it more general for any channel numbers
%  set the parameters in MAR structure from coeffient matrix A, C
%
% Usage:
%  [mar] = mar_init(A, C);
% 
%
% struct MAR:
%     mar.order               order of model
%     mar.lag(k).a            coeff matrix at lag k
%     mar.noise_cov           noise covariance matrix
%

[row, col] = size(A);  % row = chans

% set the parameters in MAR structure
mar.order = col/row;
mar.noise_cov = C;

for k=1:mar.order,
  mar.lag(k).a = A(:,1+(k-1)*row:k*row);
end




