function r = wz_boundExponential(mu, sizeOut, r1, r2);
% see: http://stackoverflow.com/questions/12358860/exponential-random-numbers-with-a-bound-in-matlab
%
% This function realizes a exponential distribution within given upper and lower boundaries.

if(~exist('r1','var') || isempty(r1))
    r1 = 0;
end

if(~exist('r2','var') || isempty(r2))
    r2 = inf;
end

minE = exp(-r1/mu);
maxE = exp(-r2/mu);

randBounded = minE + (maxE-minE).*rand(sizeOut);
r = -mu .* log(randBounded);