function chartview(dat,channel,triali,trialj,pointi,pointj)
% CHARTVIEW Chart view of data set
% 
% Input(s):
%   dat: data set in Matlab format
%   triali: start trial
%   trialj: ending trial
%   pointi: start point
%   pointj: ending point
% 
% Example:
%   chartview(data,9,1,5,1,15)
% 
% See also: sigplot, gridview.

% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.2$ $Date: 12-Sep-2007 16:41:40$
% SHIS UT-Houston, Houston, TX 77030, USA.
% 
% Lei Xu, Hualou Liang

tr1 = triali;
tr2 = trialj;
po1 = pointi;
po2 = pointj;

figure('Name','Chart View','NumberTitle','off')
for i = tr1:tr2
    subplot(tr2-tr1+1,1,i-tr1+1);
    h = gca;
    x = dat(po1:po2,channel,i);
    plot(po1:po2,x);
    set(h,'XLim',[po1 po2]);
    set(h,'XTick',po1:po2)
    ylabel(num2str(i));
    if i == tr1
        title(sprintf('Channel %d',channel));
    end%if
    if i == tr2
        xlabel('Time points');
    end
end

end%funciton

% [EOF]
    