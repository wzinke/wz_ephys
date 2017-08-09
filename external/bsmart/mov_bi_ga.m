function [Fxy,Fyx] = mov_bi_ga(dat,startp,endp,win,order,fs,freq)
% MOV_BI_GA Compute the granger causality from the moving window Bivariate models
% 
% Usage:
%   [Fxy,Fyx] = mov_bi_ga(dat,startp,endp,window,order,fs,freq) 
% 
% Input(s):
%   dat     - data set in Matlab format
%   starp   - start position
%   endp    - ending position
%   win     - window length
%   order   - model order
%   fs      - Sampling rate
%   freq    - a vector of frequencies of interest, usually freq = 0:fs/2
% 
% Output(s):
%   Fx2y    - the causality measure from x to y
%   Fy2x    - causality from y to x
%             The order of Fx2y/Fy2x is 1 to 2:L, 2 to 3:L,....,L-1 to L, where
%             L is the number of channels. That is, 1st column: 1&2; 2nd:
%             1&3; ...; (L-1)th: 1&L; ...; (L(L-1))th:(L-1)&L.
% 
% Example:
%   [Fxy,Fyx] = mov_bi_ga(data,1,18,10,5,200,[1:100])
% 
% See also: one_bi_ga.

% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.2$ $Date: 12-Sep-2007 15:53:00$
% SHIS UT-Houston, Houston, TX 77030, USA.
% 
% Lei Xu, Hualou Liang

trail   = size(dat,3);
channel = size(dat,2);
points  = endp-startp+1;
   
fxy = [];
fyx = [];
b = zeros(2,trail*win);
count = 0;
total = (points-win+1)*channel*(channel-1)/2;
hw = waitbar(0,'Please wait...');
set(hw,'Name','Granger Causality Calculation');
for ii = 1:points-win+1
    a = dat(ii+startp-1:ii+startp+win-2,:,:);
    for i = 1:(channel-1)
        for j = (i+1):channel
            waitbar(count/total);
            for k = 1:trail
                b(1,(k-1)*win+1:k*win) = a(:,i,k);
                b(2,(k-1)*win+1:k*win) = a(:,j,k);
            end
            i1 = i;
            i2 = i;
            kk = 0;
            m = channel-1;
            while i1 > 1
                kk = kk+m;
                m = m-1;
                i1 = i1-1;
            end
            kk = kk+j-i2;
            [pp,cohe,fxy(kk,:,ii),fyx(kk,:,ii)] = pwcausal(b,trail,win,order,fs,freq); %#ok<ASGLU,AGROW>
            count = count+1;
        end
    end
end
close(hw);

Fxy = fxy;
Fyx = fyx;

end%mov_bi_ga

% [EOF]