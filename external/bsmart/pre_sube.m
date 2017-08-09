function data = pre_sube(dat)
% PRE_SUBE Subtract the ensemble mean
% 
% Usage:
%  data = pre_sube(dat);
% 
%  dat: input file;
%  data: output file after subtract the ensemble mean

% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.2$ $Date: 14-Sep-2007 15:35:25$
% SHIS UT-Houston, Houston, TX 77030, USA.
%

[points channel trail] = size(dat);
for i = 1:channel
    b = dat(:,i,:);
    c = mean(b,3);      % channel ensemble mean
    for j = 1:trail
        dat(:,i,j) = dat(:,i,j)-c;  % subtract esnemble mean
    end
end
data = dat;

end%pre_sube

% [EOF]
