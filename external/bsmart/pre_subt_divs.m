function [data] = pre_subt_divs(dat)
% PRE_SUBT_DIVS Subtract the temporal mean and divide 
%               by standard deviation
% 
% Usage: 
%  [data] = pre_subt_divs(dat);
% 
%  dat: input file
% 
%  data: output file after subtract the temporal mean 
%        and divide by standard deviation

% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.1$ $Date: 11-Sep-2007 22:09:38$
% SHIS UT-Houston, Houston, TX 77030, USA.

points = size(dat,1);
channel = size(dat,2);

for i = 1:channel
    b = dat(:,i,:);
    c = mean(b,1);
    e = std(b,0,1);
    for j = 1:points
        dat(j,i,:) = dat(j,i,:)-c;
        dat(j,i,:) = dat(j,i,:)./e;
    end
end

data = dat;

end%pre_subt_divs

% [EOF]