function gridview(dat,trial,chs,che,pts,pte)
% GRIDVIEW Grid view of data set
%
%   dat: data set in Matlab format
%   trial: specified trial
%   chs: start channel
%   che: ending channel
%   pts: start point
%   pte: ending point
%
% Example:
%   gridview (data,1,6,5,1,15)

% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.2$ $Date: 12-Sep-2007 16:41:40$
% SHIS UT-Houston, Houston, TX 77030, USA.
%
% Lei Xu, Hualou Liang

% Grid view one
figure('Name','Grid View 1','NumberTitle','off')
x = dat(pts:pte,chs:che,trial);
imagesc(pts:pte,chs:che,x');
xlabel('Time points');
ylabel('Channel indexes');
tstr = sprintf('Trial %d',trial);
title(tstr);
colorbar;

% Grid view two
figure('Name','Grid View 2','NumberTitle','off');
for i = chs:che
    subplot(che-chs+1,1,i-chs+1);
    x = dat(pts:pte,i,trial);
    plot(pts:pte,x);
    axis tight;
    ylabel(num2str(i));
    if i == chs
        title(tstr);
    end%if
    if i == che
        xlabel('Time points');
    end
end

end%gridview

% [EOF]