function [Y] = mar_gen (mar, T, discardnum)
%
%   [Y] = mar_gen (mar, T)
%  Generate T samples of time series from MAR model
%   Y_t + A1*Y_(t-1) + ... =  C*W_t
% 
% Input: 
%   mar.lag(k).a   is ar coefficient matrix at lag k
%   mar.noise_cov  estimated noise covariance
%   T              number of samples to generate
%   discardnum:  discard the transient data length
%
% Output:
%      Y:  dim x T              
%

if nargin < 3,
  discardnum = 10^3;
end

p = mar.order;  % Order of AR model
d=size(mar.lag(1).a,1);
Y=zeros(d,T+discardnum);

randvec = gaussian(T+discardnum,zeros(1,d),mar.noise_cov)';  % dim x T+Ndisc
% Generate first p elements
Y(:,1:p) = randvec(:,1:p);

% Generate rest of series
for i=p+1:T+discardnum,
  for k=1:p,
    Y(:,i) = Y(:,i) - mar.lag(k).a*Y(:,i-k);
  end
  % Y(:,i)=Y(:,i)+gaussian(1,zeros(1,d),mar.noise_cov)';
  Y(:,i) = Y(:,i) + randvec(:,i);
end

Y = Y(:,discardnum+1:discardnum+T);








