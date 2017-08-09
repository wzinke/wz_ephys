function ratio = consistencytest(dat,arcoeff,arnoise)
% CONSISTENCYTEST Consistency test
% 
% Sytax:
%   ratio = consistencytest(data,A,Ve)
% 
%   dat: data set in Matlab format
%   arcoeff: AR coefficient 
%   arnoise: AR noise
% 
% Ouput(s): 
%   ratio in Matlab format
% 
% Example: 
%   

% Copyright (c) 2006-2007 BSMART Group
% by Richard Cui
% $Revision: 0.2$ $Date: 15-Sep-2007 09:04:13$
% SHIS UT-Houston, Houston, TX 77030, USA.
% 
% Lei Xu, Hualou Liang

points = size(dat,1);
trial = size(dat,3);

aredat = arcoeff;
arndat = arnoise;
s1 = size(arndat);
channel = sqrt(s1(2));
s2 = size(aredat);
window = points+1-s2(1);

[R1,R2,diffR,ratio] = movingwin(dat,aredat,arndat,window,points,trial,channel);

ratio = 100*(1-ratio);

end%function

% [EOF]