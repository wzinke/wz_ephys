function conetwork(coherence,location,thre,time,fre1,fre2,chan)
% CONETWORK Analysis and visualization of the coherence network
% 
%   coherence: pair coherence
%   location: location of the sites
%   thre: threshold
%   time: specify the window number
%   fre1: starting frequency bin
%   fre2: ending frequency bin
%   chan    - Channel of interest
% 
% Example:
%   conetwork(coherence,location,0.25,5,1,50,[9,10,13]);
% 
% See also: ganetwork.


% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.4$ $Date: 15-Sep-2007 19:12:50$
% SHIS UT-Houston, Houston, TX 77030, USA.
%
% Lei Xu, Hualou Liang

% parameter settings
t   = time;
x1  = fre1;
x2  = fre2;
dat = coherence;
dat2= location;
s1  = size(dat);

% error checking
if (x1 > s1(2))
    errordlg('please input correct start frequency','parameter lost');
    return
end
if (x2 > s1(2))
    errordlg('please input correct end frequency','parameter lost');
    return
end
if (t > s1(1))
    errordlg('please input correct window','parameter lost');
    return
end

% find the channel number
N = s1(3);
channel = (1+sqrt(1+8*N))/2;    % k = (1+sqrt(1+8N))/2, N = number of pairs of coherence

% create channel labels
label = cell(1,channel);
for i = 1:channel
    label(1,i)={num2str(i)};
end
% control circle type
circle = zeros(channel,1);  

% dag & dag2 for the function 'drawlayout'
dag = zeros(channel,channel);
dag2= dag;
left = dat2(1,:);
right = dat2(2,:);
data = squeeze(dat(t,:,:));
xlabel = x1:x2;
en = channel*(channel-1)/2;
for i = 1:en
    peak  = fpeak(xlabel,data(xlabel,i),30,[x1,x2,0,1]);
    sizep = size(peak);
    if sizep(1)>=1
        for pi=1:sizep(1)
            if peak(pi,2)>thre
                kk=i;
                m=channel-1;
                while kk>0
                    kk=kk-m;
                    m=m-1;
                end
                ii=channel-1-m;
                jj=ii+kk+m+1;
                dag(ii,jj)=1;
                dag(jj,ii)=1;
            end
        end
    end
end
% set channels of interets
tp = true(1,channel);
tp(chan) = false;
chanary = 1:channel;
notchan = chanary(tp);
circle(chan) = 1; 
dag(notchan,:) = 0;
dag(:,notchan) = 0;
% plot figure
figure('Name','Coherence Network','NumberTitle','off')
draw_layout(dag,label,circle,left,right,'ch');
tstr = sprintf('Window: %d, Frequency: %.0f - %.0f Bins, Threshold: %.2f',time,fre1,fre2,thre);
title(tstr);

% if no peak is found
if dag==dag2
    for i=1:en
        for j=x1:x2
            if data(j,i)>thre
                kk=i;
                m=channel-1;
                while kk>0
                    kk=kk-m;
                    m=m-1;
                end
                ii=channel-1-m;
                jj=ii+kk+m+1;
                dag(ii,jj)=1;
                dag(jj,ii)=1;
            end
        end
    end
    % set channels of interets
    tp = true(1,channel);
    tp(chan) = false;
    chanary = 1:channel;
    notchan = chanary(tp);
    circle(chan) = 1;
    dag(notchan,:) = 0;
    dag(:,notchan) = 0;
    % plot figure
    figure('Name','Coherence network','NumberTitle','off')
    draw_layout(dag,label,circle,left,right,'ch');
    tstr = sprintf('No peak found. Window: %d, Frequency: %.1f - %.1f Hz, Threshold: %.2f',time,fre1,fre2,thre);
    title(tstr);
end

end%function

% [EOF]