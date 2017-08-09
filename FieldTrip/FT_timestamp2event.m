function [plxmat, cfg] = FT_timestamp2event(cfg, evtimes, duration)
% FT_TIMESTAMP2EVENT takes a vector of time stamps and converts them into
% an event structur array that will be appended to the Fieldtrip cfg struct.
%
% In Fieldtrip vents are represented as
%   event.type      string
%   event.sample    expressed in samples, the first sample of a recording is 1
%   event.value     number or string
%   event.offset    expressed in samples
%   event.duration  expressed in samples
%   event.timestamp expressed in timestamp units, which vary over systems (optional)
%
%
% wolf zinke, 22.4.2016

evtimes(isnan(evtimes)) = []; % remove NAN entries

numEv = length(evtimes);

if(numEv == 0)
    % nothing to do
    disp('No event times to convert!');
    return;
end

if(~exist('duration','var') || isempty(duration))
    duration = repmat(1,numEv,1);
end


% The sample numbers returned in event.sample correspond with the
% timestamps, correcting for the difference in sampling frequency in the
% continuous LFP channels and the system sampling frequency. Assuming 40kHz
% sampling frequency for the system and 1kHz for the LFP channels, it is
%   event.sample = timestamp / (40000/1000);
 event.sample = timestamp / (40000/1000);