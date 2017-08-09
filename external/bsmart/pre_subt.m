function data = pre_subt(dat)
% PRE_SUBT Subtract the temporal mean
% 
% Usage: 
%  data = pre_subt(dat);
% 
%  dat: input file;
%  data: output file after subtract the temporal mean

% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.1$ $Date: 11-Sep-2007 22:09:38$
% SHIS UT-Houston, Houston, TX 77030, USA.
% 
% Lei Xu, Hualou Liang

points = size(dat,1);
channel = size(dat,2);

for i = 1:channel
    b = dat(:,i,:);
    c = mean(b,1);
    for j = 1:points
        dat(j,i,:) = dat(j,i,:)-c;
    end
end

data = dat;

end%pre_subt

% [EOF]