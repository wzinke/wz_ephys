function ns = NeuroSound(SS,tbound,pbspeed,ampscale,basefreq)
% NEUROSOUND audio output event-based recording that really makes a
% SqueakySpk Squeak.
%
% NS = NEUROSOUND(SS,PBSPEED,AMPSCALE,BASEFREQ)
%
%       Inputs:
%       SS = SqueakySpk object
%       TBOUND = [t1 t2] defines the time window, in seconds, of the
%       recording that you want to hear.
%       PBSPEED = scale factor for the playback speed. May be helpful in
%       detecting patterns that occur on different time scales. Increase to
%       make the playback faster. PBSPEED = 1 is real-time.
%       AMPSCALE = Scale factor to keep sound signal within the dynamic
%       range of your speakers [+-1 volt]. Adjust up to stretch dynamic
%       range over [+-1 volt] (e.g. you can't discriminate events in the
%       recording from background noise) and down if you are getting
%       saturation. This option is available only for
%       continuous waveforms and will be ignored otherwise.
%       BASEFREQ = Lowest frequency to make make a scale of notes
%       corresponding to units or channels
%
%       Outputs:
%       NS = Neurosound audio object. This oject has the following methods
%       associated with it:
%           play(ns) - play the neurosound
%           stop(ns) - stop the neurosound
%           pause(ns) - pause the neurosound
%           resume(ns) - resume the neurosound after a pause
%
%       Created by: Jon Newman (jnewman6 at gatech dot edu)
%       Location: The Georgia Institute of Technology
%       Created on: Apr 24, 2010
%       Last modified: Apr 24, 2010
%
%       Licensed under the GPL: http://www.gnu.org/licenses/gpl.txt

% check number and type of arguments
if( nargin < 3)|| isempty(tbound)
    tbound = [];
end
if( nargin < 3)|| isempty(pbspeed)
    pbspeed = 1;
end
if(nargin < 4) || isempty(ampscale)
    ampscale = 1;
end
if(nargin < 5) || isempty(basefreq)
    basefreq = 110;
end

if isempty(SS.unit)
    chan = SS.channel(SS.clean);
    spktime = SS.time(SS.clean);
    tonerep = 'channels';
else
    chan = SS.unit(SS.clean&SS.unit~=0);
    spktime = SS.time(SS.clean);
    tonerep = 'units';
end

% Create continuous waveform from the discrete events
if isempty(tbound)
    sortedtimes = sort(spktime); % sorted times, in seconds, referenced to t=0.
else
    sortedtimes = sort(spktime(spktime>tbound(1)&spktime<tbound(2)));% sorted times, in seconds, referenced to t=0.
    sortedtimes = sortedtimes - min(sortedtimes);
end
sortedtimes(sortedtimes==0) = []; % get rid of any time that is zero
maxT = sortedtimes(end) + 1; % length of sound will be the time of the final event + 1 sec

N = ceil(SS.fs*maxT);
wave = zeros(N,1);

eventInd = ceil(SS.fs*sortedtimes);

% Create 10 ms tone snips for each channel
numchan = max(chan);
snip = zeros(numchan,length(1/SS.fs:1/SS.fs:pbspeed*0.1));
for k = 1:numchan
    snip(k,:) = sin((basefreq*2^(k-1)/12)/pbspeed*2*pi*(1/SS.fs:1/SS.fs:pbspeed*0.1));
end

for k = 1:length(eventInd)
    wave(eventInd(k):eventInd(k)+length(snip)-1) = snip(chan(k),:);
end

ns = audioplayer(wave,SS.fs*pbspeed);

% Instructions for use
disp(' ');
disp('You have created a neurosound object with the following properties:');
disp(['Playback speed: ' num2str(pbspeed) ' X real-time']);
disp(['Amplitude scale factor: ' num2str(ampscale)]);
disp(['Tones represent different' tonerep])
disp(' ');
disp('Use the following commands to control playback: ');
disp('[1]: "play(ns)" to stop the playback');
disp('[2]: "stop(ns)" to stop the playback');
disp('[3]: "pause(ns)" to pause the playback');
disp('[4]: "resume(ns)" to resume the playback');

end

