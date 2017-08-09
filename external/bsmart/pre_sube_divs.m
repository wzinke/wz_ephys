function [data] = pre_sube_divs(dat)
% PRE_SUBE_DIVS Subtract the ensemble mean and divide by standard deviation
%
% Usage:
%  [data] = pre_sube_divs(dat);
%
%  dat: input file;
%  data: output file after subtract the ensemble mean and divide by
%  standard deviation

% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.1$ $Date: 11-Sep-2007 22:29:45$
% SHIS UT-Houston, Houston, TX 77030, USA.


[points channel trail] = size(dat);
for i = 1:channel
    b = dat(:,i,:);
    c = mean(b,3);
    e = std(b,0,3);
    for j = 1:trail
        dat(:,i,j) = dat(:,i,j)-c;
        dat(:,i,j) = dat(:,i,j)./e;
    end
end
data=dat;

end%pre_sube_div

% [EOF]