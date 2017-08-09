function po_view(power,fs,chx,timen)
%  PO_VIEW View power
% 
%   power   - auto power data set
%   fs      - sampling rate
%   chx     - specified channel
%   timen   - (optional) view coherence at one time
% 
%  Example:
%   po_view(autospect,200,9); 
%   po_view(autospect,200,9,2);
% 
% See also: autopower, bi_power.

% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.2$ $Date: 12-Sep-2007 15:53:00$
% SHIS UT-Houston, Houston, TX 77030, USA.
%
% Lei Xu, Hualou Liang

dat = power;
if nargin < 4
    si = size(dat);
    tstr = sprintf('Autopower of Channel %d (dB)',chx);
    if si(1) == 1   % only one window
        c = dat(1,:,chx);
        c = 10*log10(c);    % in dB
        freq = linspace(0,fs/2,si(2));
        figure;
        plot(freq,c);
        h = gca;
        title(tstr);
        xlabel(h,'Frequency (Hz)')
        ylabel(h,'Auto Power (dB)')
    else % more than one window
        b = dat(:,:,chx);
        b = 10*log10(b);
        time = [1 si(1)];
        freq = [0 fs/2];
        figure;
        imagesc(time,freq,b');
        colorbar;
        axis xy;
        h = gca;
        title(tstr);
        xlabel(h,'Time')
        ylabel(h,'Frequency (Hz)')
    end
else % window specified
    si = size(dat);
    c = dat(timen,:,chx);
    c = 10*log10(c);
    freq = linspace(0,fs/2,si(2));
    figure;
    plot(freq,c);
    h = gca;
    tstr = sprintf('Autopower of Channel %d at Window %d',chx,timen);
    title(tstr);
    xlabel(h,'Frequency (Hz)')
    ylabel(h,'Autopower (dB)')
end

end%function

% [EOF]