function [W, X, par, parName] = waveformPreprocessing(waveforms,cfg,isolationquality)
% function [W, X, par, parName] = waveformPreprocessing(cfg,waveforms,isolationquality)
% --- aligns, interpolates, normalizes waveforms of units and measures their waveform parameters parName
% using isolationquality is optional

% default behavior
if ~isfield(cfg, 'troughalign'),    cfg.troughalign = 1;        end % alignment at waveform trough
if ~isfield(cfg, 'interpn'),        cfg.interpn = 10;           end % x10
if ~isfield(cfg, 'interptype'),     cfg.interptype = 'splines'; end   % type of interp: 'splines', 'pchip', 'linear' % used splines by default
if ~isfield(cfg, 'normalize'),      cfg.normalize = 1;          end % between -1 and 1
if ~isfield(cfg, 'fsample'),        cfg.fsample = 40000;        end % 40 kHz

[numUnits numSamples] = size(waveforms);

W = waveforms;
X = 0:numSamples-1;

% --- interpolation
if cfg.interpn>1
    Xn = linspace(X(1), X(end), length(X)*cfg.interpn);
    Wn= zeros(numUnits,length(Xn));
    for k=1:numUnits
        if(~all(isnan(W(k,:))))
            indnan=isnan(W(k,:));
            tmpW=W(k,:);
            tmpW=tmpW(~isnan(tmpW));
            tmpX=X(~isnan(tmpW));
            Wn(k,:) = interp1(tmpX,tmpW,Xn,cfg.interptype,nan);
        end
    end,
    W = Wn;
    X = Xn;
end

if cfg.troughalign
    % --- make new xscale (as short as possible)
    [~,b]=min(W,[],2);
    lenL = 1-max(b);
    lenR = length(X)-min(b);
    Xn = lenL:lenR;
    % --- shift waveforms trough to zero
    Wn=nan(numUnits,length(Xn));
    for k=1:numUnits
        Wn(k,(max(b)-b(k))+(1:length(X)))=W(k,:);
    end
    W = Wn;
    X = Xn;
end

% --- X converted from interpolated samples to time (in ms)
X = X/(cfg.fsample*0.001*cfg.interpn); % time steps of 2.5 microseconds: 1/(40kHz * 0.001 in ms * 10 from interp)

% --- normalize waveforms between -1 and 1
if cfg.normalize,
    for k=1:numUnits
        if(max(W(k,:))~=min(W(k,:)))
            W(k,:) = (W(k,:) - min (W(k,:)))/(max(W(k,:)) - min (W(k,:)));
            W(k,:) = 2*W(k,:)-1;
%              if (k==1)
%                  figure
%                  plot(X,W(k,:),'.')
%                  return
%              end
        end
    end
end

% --- extract waveform parameters: peak to trough and time for (25%) repolarization

par = nan(numUnits,3);

for k=1:numUnits

    Wk = W(k,:);

    if ~any(~isnan(Wk))
        continue
    end

    if nargin==3 % using isolationquality is optional
        if isolationquality(k)~=3
            continue
        end
    end

    % peak to trough duration

    % best method to get the right values
    iminima=minima(Wk,1); % getting all minima, it's assured that 1 element around isn't nan
    [~,iTmp]=min(Wk(iminima));
    idown=iminima(iTmp);
    [~,iup]=max(Wk(idown+1:end));

    if isempty(iup) || isempty(idown)
        continue
    end

    iup=idown+iup;

    % maximum rate of fall/maximum rate of rise
    % differentiate first:

    wdiff = diff(Wk,1); %last ind is original ind-1
    %find the global minimum before the action potential through
    maxrise = min(wdiff(1:idown));
    %find the global maximum after the action potential through
    maxfall = max(wdiff(idown+1:iup));

    risefallratio = maxfall/maxrise;

    par(k,1)=X(iup)-X(idown);
    if(par(k,1)<0) % checks to avoid bad surprises
        display('ERROR: negative peak-to-trough duration in the trial')
        display(k)
    end
    parName{1}='Peak to trough duration';

    % amplitude of the peak before depolarization (checked but not used)

    tmp = find(~isnan(Wk));

    if(~isempty(iminima(iminima<idown)))
        [~,iTmpMin]=min(Wk(iminima<idown));
        idownPrePeak=iminima(iTmpMin);
    else
        idownPrePeak=tmp(1);
    end

    % getting all maxima before the minimum
    imaxima=maxima(Wk(1:idown),1); % it's assured that 1 element around isn't nan
    if(isempty(imaxima))
        par(k,2)=0;
    else
        [~,iTmpMax]=max(Wk(imaxima));
        iupPrePeak=imaxima(iTmpMax);
        if(iupPrePeak<=tmp(1))
            par(k,2)=0;
        else
            par(k,2)=Wk(iupPrePeak)-Wk(idownPrePeak);
            if(par(k,2)<0) % checks to avoid bad surprises
                display('ERROR: negative amplitude for the pre-peak in trial')
                display(k)
            end
        end
    end
    parName{2}='Amplitude of peak before depolarization';

    % time for repolarization (25% of peak to trough)
    percentage=25;

    amplPeaktoTrough=Wk(iup)-Wk(idown);
    repInd=nearest(Wk(iup+1:end),Wk(iup)-percentage/100*amplPeaktoTrough); % index of the nearest value to the decay
    if(~isempty(repInd)) % checking it is not empty
        if(~isnan(Wk(iup+repInd+1))) % checking it is not the last value of the serie before all nans
            par(k,3)=X(iup+repInd)-X(iup);
        else
            par(k,3)=nan;
        end
    end
    parName{3}='Time for repolarization';
end
