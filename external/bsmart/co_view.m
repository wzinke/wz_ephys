function co_view(coh,chx,chy,fs,timen)
% CO_VIEW  View coherence
%
% Syntax:
%   co_view(coh,chx,chy,fs)
%   co_view(coh,chx,chy,fs,timen)
%
% Input(s):
%   coh     - Coherence data set
%   fs      - Sampling rate (Hz)
%   chx     - One channel
%   chy     - The other channnel
%   timen   - A 2D array of window index (optional)
%             [index1,center_time1;index2,center_time2;...;indexk,center_timek]
%
% Example:
%   co_view(paircoh,9,10,200);
%   co_view(paircoh,9,10,200,[1,30]);

% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.3$ $Date: 13-Sep-2007 10:08:16$
% SHIS UT-Houston, Houston, TX 77030, USA.
%
% Lei Xu, Hualou Liang

% parameter settings
[nw,nb,c] = size(coh);  % nw = number of windows, nb = No. frequency bins, c = coh pairs #
N = (1+sqrt(1+8*c))/2;  % total number of channels
if chx == chy
    error('Coherrence must be considered between different channels.');
elseif chx > chy        % swap them
    t = chx;
    chx = chy;
    chy = t;
end%if
k = (2*N-chx)*(chx-1)/2+(chy-chx);  % note: chx < chy
cohk = squeeze(coh(:,:,k));         % coherence for specific channel pair

% plot
if nw == 1
    if nargin <5
        % plot spectrum, not specified
        timen = [];
    else % nargin == 5
        if size(timen,1) > 1    % more than one window
            error('Cannot specify more than one window');
        end%if
    end%if
    drawspect(cohk,fs,timen);  % draw spectrm
else % nw >= 2
    if nargin == 5
        % plot grid view and spectra
        drawspect(cohk,fs,timen);
        time = timen(:,2);
        if length(time) == nw
            drawgrid(cohk,fs,time,chx,chy);
        end%if
    elseif nargin < 5
        % plot grid view only
        time = 1:nw;
        drawgrid(cohk,fs,time,chx,chy);
   end%if
end%if

end%function

% -------------
% subroutines
% -------------
function drawspect(coh,fs,timen)

% coh = windows x frequency

[nw,nbin] = size(coh);

if isempty(timen)
    winc = 1:nw;
else
    winc = 1:size(timen,1);
end%if
f = linspace(0,fs/2,nbin);

% draw spectrum window by window

for k = winc
    if isempty(timen)
        spect = coh;
    else
        spect = coh(timen(k),:);
    end%if
    figure
    plot(f,spect)
    xlabel('Frequency (Hz)')
    ylabel('Coherence')
    if ~isempty(timen)
        title(sprintf('Window centered at %d ms',timen(k,2)));
    end%if
end%for

end%function

% -------------
function drawgrid(coh,fs,time,chx,chy)

[nw,nbin] = size(coh);
f = linspace(0,fs/2,nbin);

figure
imagesc(time,f,coh')
axis xy
colorbar

xlabel('Time')
ylabel('Frequency (Hz)')
tstr = sprintf('Coherence between Channel %d and Channel %d',chx,chy);
title(tstr);

end%function

% [EOF]