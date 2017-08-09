function [out_val, parvec] = wz_gamma_mod(X, amp, kapa, tau, lag, offset)
% gamma distribution.
%
% wolf zinke, 24.10.2004 (modified May 2014)

npar = nargin -1;
% input comes in the form of (Xvals, [coeffarr])
if(npar == 1)
    parvec = amp;
    npar = length(parvec);
    amp  = parvec(1);
    kapa = parvec(2);
    tau  = parvec(3);
    if(np > 3)
        lag = parvec(4);
    else
        lag = 0;
    end
    if(np > 4)
        offset = parvec(5);
    else
        offset = 0;
    end
end

if(~exist('lag','var') || isempty(lag))
    lag = 0;
end

if(~exist('offset','var') || isempty(offset))
    offset = 0;
end

if(lag ~= 0)
    pos = find(abs(X-lag) == min(abs(X-lag)));
    tmp_out = zeros(1,length(X));
    X = X(1:end-pos);
end

out_val = X .^(kapa-1).*exp(-X./tau)./gamma(kapa)./(tau.^kapa);

out_val = amp * (out_val/max(out_val)) + offset;

if(lag ~= 0)
    tmp_out(pos+1:end) = out_val;
    out_val = tmp_out;
end

if(npar == 5)
    parvec = [amp, kapa, tau, lag, offset];
elseif(npar == 4)
    parvec = [amp, kapa, tau, lag];
else
    parvec = [amp, kapa, tau];
end
    

