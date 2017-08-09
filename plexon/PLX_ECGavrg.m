function [peakpos] = PLX_ECGavrg(cchn, lbl, peakpos)
wvwin = [-650, 650];

cchn(isnan(cchn)) = [];

if(~exist('peakpos','var') || isempty(peakpos))
    % % apply a high pass filter
    hpFilt = designfilt('highpassiir','FilterOrder',20, ...
                        'PassbandFrequency',0.5,'PassbandRipple',0.1, ...
                        'SampleRate',1000);
   cchnfilt = filter(hpFilt,cchn);

y = sgolayfilt(cchnfilt, 0, 21);

    [peakpos] = peakdet(y, mad(y,1)/2);
else
    y = cchn;
end

ecgwave = nan(length(peakpos),diff(wvwin)+1);

for(c=1:length(peakpos))
    if(peakpos(c) > abs(wvwin(1)) && peakpos(c) < length(y)-wvwin(2))
        ecgwave(c,:) = y(peakpos(c)+wvwin(1):peakpos(c)+wvwin(2));
    end
end

meanwave = nanmean(ecgwave);

Yrng = [min(meanwave), max(meanwave)];
plot([0 0], Yrng, 'r');

plot(wvwin(1):wvwin(2), meanwave);

ylim(Yrng);
xlim(wvwin);

if(exist('lbl','var') && ~isempty(lbl))
    text(0.975*wvwin(2), Yrng(1),  lbl, 'FontSize', 12, 'Interpreter', 'none', ...
    'HorizontalAlignment','right','VerticalAlignment','bottom', 'Rotation',0);
end