function LE = lyap_batch(arcoeff,arnoise,T)
% LYAP_BATCH A batch version to compute Lyapunov exponen
% 
% Usage: 
%   LE = lyap_batch(arcoeff,arnoise,T);
% 
%   T   -   Number of samples to generate;
% 
% Note: 
%   In fact we do not use ARNOISE in Lyapunov calculation, 
%   but need it in MAR_make routine
%

% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.2$ $Date: 12-Sep-2007 10:39:05$
% SHIS UT-Houston, Houston, TX 77030, USA.
%
% Hualou Liang, 02/09/99, FAU
% 

LE = [];
len = size(arnoise,1);
h = waitbar(0,'Please wait...');    % TODO: set a window title
for t = 1:len,
    waitbar(t/len);
    mar = MAR_make(arcoeff(t,:), arnoise(t,:));
    le = lyap(mar,T,1000);          % LE = lyap(mar, T, discardnum);
    LE = cat(2,LE,le);
end;
close(h);

end%lyap_batch

% [EOF]
