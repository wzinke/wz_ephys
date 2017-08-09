function [LE] = lyap (mar, T, discardnum)

% calculating Lyapunov Exponent of MVAR
% if the exponent > 0. the process is not stationary
% Using L1, L2 distance measure, result are same.
%   Y_t + A1*Y_(t-1) + ... =  0
%  Generate T samples of time series from MAR model
%   Y_t + A1*Y_(t-1) + ... =  C*W_t
% 
% Input: 
%   mar.lag(k).a   is ar coefficient matrix at lag k
%   mar.noise_cov  estimated noise covariance(do not need it)
%   T              number of samples to generate
%   discardnum:  discard the transient data length
%
% Output:
%      Le:  Lyapunov exponent, a scalar number
%

% 
% Hualou Liang, 02/08/99, FAU
%



LE = 0;
if nargin < 3,
  discardnum = 10^3;
end

p = mar.order;  % Order of AR model
d=size(mar.lag(1).a,1);
Y=zeros(d,T+discardnum);

randvec = rand(d, p);  % dim x p 
% Generate first p elements
% Y(:,1:p) = randvec(:,1:p);

% Y(:,1:p) = randvec ./ (ones(d,1)*sqrt(sum(randvec.^2)));  % L2
% Y(:,1:p) = randvec ./ (ones(d,1)*sum(abs(randvec)));     % L1

%%%% change HERE  %%%%%%%%%%%%
Y(:,1:p) = randvec / sum(sum(abs(randvec)));     % L1
% Y(:,1:p) = randvec / sqrt(sum(sum(randvec.^2)));     % L2

% Generate rest of series to be discarded
for i=p+1:discardnum,
  for k=1:p,
  	Y(:,i) = Y(:,i) - mar.lag(k).a*Y(:,i-k);
  end
  
  % d = sum(abs(Y(:,i)));
  % d = sqrt(sum(Y(:,i).^2));
  % Y(:, i) = Y(:,i)/d;

%%%% change HERE  %%%%%%%%%%%%
  Y(:,i-p+1:i) = Y(:,i-p+1:i) / sum(sum(abs(Y(:,i-p+1:i))));   % L1
%  Y(:,i-p+1:i) = Y(:,i-p+1:i) / sqrt(sum(sum(Y(:,i-p+1:i).^2)));  % L2
  
end

% for i=p+1:T+discardnum,
%d0 = sqrt(sum(Y(:,discardnum).^2));

for i=discardnum+1:T+discardnum,
  for k=1:p,
  	Y(:,i) = Y(:,i) - mar.lag(k).a*Y(:,i-k);
  end
  
  % d = sum(abs(Y(:,i)));       % L1
  % d = sqrt(sum(Y(:,i).^2));   % L2
  % Y(:, i) = Y(:,i)/d;

%%%% change HERE  %%%%%%%%%%%%
  d = sum(sum(abs(Y(:,i-p+1:i))));
%  d = sqrt(sum(sum(Y(:,i-p+1:i).^2))); 
  Y(:,i-p+1:i) = Y(:,i-p+1:i) / d;

  LE = LE + log(d);
  % pause
end

LE = LE/T;

% d1 = sqrt(sum(Y(:,discardnum+T).^2));
% LE = log(d1/d0)/T;





