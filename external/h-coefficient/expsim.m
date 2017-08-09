function [peristim,peristim_cntr,grndtruth,grndtruth_cntr,rawdata] = expsim(duration,nstims,percent,noisefreq,noiseamp,xdata,sigma,rspdealy,amp,bslrange,mfr,doplot)



% *********************************************************************************************
% *********************************************************************************************
% *******   This code was authored by Michael Hill (buffalohill@me.com), it may not     *******   
% *******   be used as a whole or in parts without the written permission of Michael    *******  
% *******   Hill and any use of this code (in whole or in parts) must reference the     ******* 
% *******   paper in which the h-coefficient method was first introduced (contact me    ******* 
% *******   to inquire about appropriate referencing formats for this code).             ****** 
% *********************************************************************************************
% *********************************************************************************************
%
%
% expsim
% expsim is a script that generates artificial spiking data for a whole
% experimental session, with sections cut out from that session, which
% may be used to plot rasterplots or peristimulus time histograms (PSTH).
%
% Input:
% expsim takes the following input parameters and their values:
%
% duration   **********************************************************************************
% expsim generates a train of spikes from one artificial neuron for the
% whole duration of your experiment with responses occuring after your
% stimuli according to your response parameters set below. The input value
% 'duration' is the total duration of your experiment in minutes. For
% example you might want to generate an experiment that lasts 50 minutes
% during which you presented 50 stimuli (i.e. an interstimulus interval of
% 60 sec). Further down you define the length of your PSTH (xdata), make
% sure that your 'duration' is long enough to contain all your individual
% peristimulus response and baseline durations without overlap, i.e. your
% 'duration' must be longer then xdata * nstims (see below). 'duration' is
% entered as a scalar in units of minutes, e.g. duration = 50 (minutes).
%
% nstims   ************************************************************************************
% nstims is the number of stimuli that you present during your whole
% experiment (i.e. your number of trials). This will also reflected in the
% number of lines in your rasterplot and the confidence inrervalls in your
% PSTH. nstims is entered as a scalar in units, e.g. nstims = 50.
%
% percent   ***********************************************************************************
% A neuron might not always show the same response to the each stimulus, 
% this might be because the stimuli are different, or are interpreted
% differently or because the neuron might respond differently due to 
% implicit known or unknown parameters. 'percent' denotes the sensitivity
% of the artificial neuron to the stimuli by defining the percentage of 
% trials during which a response took place. For example:
% percent = 1       the neuron responded in each trial
% percent = 0.7     70% of the trials will show a response whilst the remainig 30% will remain flat
% percent = 0       the neuron shows no responses to the stimuli.
% percent = -0.7 	70% of the trials will show a response but in this case the remainig 30% won't be flat but will instead show a negative response of the same amplitude as amp (below) but of course never reaching below 0 Hz.
% percent is entered as a scalar, e.g. -1 =< nstims =< 1.
%
% noisefreq   *********************************************************************************
% The frequency range of a sine wave, which will be convoluted with your 
% data as a noise signal. the first value is the frequency of the noise, 
% the second value is the amount by which the frequency is randomly 
% modulated in each trial(i.e. +/- that amount).
% noisefreq is entered as a row vector with two values in units of Hz, 
% e.g. noisefreq = [5,10].
%
% noiseamp   **********************************************************************************
% The amplitude of the above noise sine wave, the first value is the 
% amplitude of the noise, the second value is the amount by which it is
% randomly modulated in each trial (i.e. +/- that amount, the actual
% values of the noiseamp are entered in multiples of the baseline value,
% see below). Be aware, expsim's main goal is to generate PSTHs for further
% analysis, the noise is added to each stimulus separatly and will cause
% some edge artifact in the raw data across the whole 'duration'. These
% artifacts can be ignored in most applications, or can even be viewed as
% additional, non-sinusoidal noise, which may make any post-hoc analysis 
% even more conservative (as is be the case in the paper published along 
% with this code). However, in other applications this noise may be 
% problematic and would have to be removed (noiseamp = [0,0]). noiseamp is 
% entered as a row vector with two values in units of Hz, e.g. noiseamp = [0.5,1].
%
% xdata   ************************************************************************************
% xdata denotes the total range of your data around each stimulus in
% milliseconds. For example, if you wish to generate PSTHs reaching from
% -5 sec to +5 sec around your stimulus, with a resolution of 1 ms (you can
% always zoom into your PSTH at a later stage, these values should be seen
% more as your maximum range) then you need to enter these values into the
% xdata variable in the format: xdata = (1:1:10000), i.e. the total xdata
% range is 10000 ms long. The script can easily be optimised to take other
% range and resolution values for xdata, but not all combinations have been
% tested and some adjustments might be necessary, especially if the
% resolution is changed. For most applications we suggest leaving xdata at
% the default values given below:  xdata = (1:1:10000). If you do change
% these values, make sure that length(xdata) is an even number so your
% stimulus can be set at zero (i.e. xdata(1,end)/2) does not produce a
% fraction.
% xdata is entered as a row vector, e.g. xdata = (1:1:10000).
%
% sigma   ************************************************************************************
% The width of the neurons response to a stimulus in ms (i.e. the standard
% deviation of your gaussian). sigma is entered as a scalar in units of
% milliseconds, e.g. sigma = 300 (ms).
%
% rspdealy   *********************************************************************************
% The location of the neuron's response peak in ms relative to the
% beginning of your xdata range. For example in the above xdata range a
% post-stimulus response with a response delay of 500 ms after stimulus
% onset would be enterred as rspdelay = 5500 (5000 ms pre-stimulus + 500
% post-stimulus time, i.e. the distance of the center of the response
% gaussian relative to the beginning of the xdata range). rspdealy is
% entered as a scalar in units of milliseconds, e.g. rspdealy = 5500 (ms).
%
% amp   *************************************************************************************
% The amplitude of your response relative to the neurons baseline firing
% rate. For example:
% amp = 0   The neuron responds with total inhibition, i.e. 0 Hz
% amp = 1   The neuron shows no response, i.e. it keeps on firing at baseline
% amp = 2 	The neuron shows a response with an amplitude double of your baseline firing rate
% amp is entered as a scalar normalized to baseline, e.g. amp = 2.
%
% bslrange   ********************************************************************************
% The beginning and end of your baseline range. This baseline range is
% the total duration before your stimulus, during which you don't expect
% any response activit in your neuron. For exmple, if you have a 5 sec
% prestimulus interval, you could enter: bslrange = [0,4800], which would
% mean that all of the 5 seconds would be used except for the last 200 ms
% where some change in firing rate might be expected if, for example, you
% set your sigma to a very large value or your rspdealy to a very low value such
% that the beginning of your response reaches beyond the stimulus itslef.
% Such settings might be appropriate in the case where a stimulus is time
% locked and therefore might be expected by the subject. bslrange needs to
% always be equal or smaller then xdata/2. bslrange should not be confused
% with baseline period. The baseline range is only used to run the Monte
% Carlo simulation the right number of times to make sure that the baseline
% firing rate of your neuron is correct. The baseline period, used
% elsewhere and mentioned in the paper, which was published along side this
% code, denotes a period, to which the response period is compared, for
% example in a t-test. Usually the baseline period will be smaller then
% bslrange and will be contained within baslrange.
% bslrange is entered as a row vector with two values in the same units as
% your xrange (usually milliseconds) e.g. bslrange = [0,4800].
%
% mfr   *************************************************************************************
% mfr is the mean firing rate during your baseline range in Hz. The Monte
% Carlo simulation for each trial will run until it has generated data with
% this baseline value over the duration of the bslrange (see above), as the
% simulated data is binary, if you want it to be precise, then make sure
% that your baseline range doesn't create fractions of Hz as it's mean
% firing rate in the following equation:
%       bsl_spikes = mfr/(1000/norm(bslrange(1,2)-bslrange(1,1)));
% where bsl_spikes will be the absolute number of spikes needed in the
% bslrange to fulfil the above mean baseline firing rate for any given
% trial.
% mfr is entered as a scalar in units of Hz, e.g. mfr = 5 (Hz).


dbstop if error

%cluster_class = [];
%rawstims = [];

%% Default settings

if nargin < 1;      
    duration    =   30;
    nstims      =   1;
    percent     =   1;
    noisefreq   =   [0,1];
    noiseamp    =   [0.5,1];
    xdata       =   (1:1:10000);
    sigma       =   300;
    rspdealy    =   5450;
    amp         =   1.25;
    bslrange    =   [0,4800];
    mfr         =   3;
    doplot      =   1;             
end

%% Create and adjust some variables

peristim = {};
rawdata = [];
duration = duration * 60 * 1000;                                            % Convert duration from minutes to ms
respstims = round(abs(percent) * nstims);                                   % Number of stimuli to which the neuron responds
bsl_spikes = mfr/(1000/norm(bslrange(1,2)-bslrange(1,1)));                  % absolute number of spikes needed in bslrange to fulfil the mean baseline firing rate (mfr).
block = duration / (nstims+2);                                              % Create the right number and duration of respons + baseline periods to fill your whole experiment duration
bslblock = block - length(xdata);
ylims = [0,1 + amp + amp*(noiseamp(1,1)+noiseamp(1,2))];

%% Create the artificail data

for a = 1:nstims                                                            % Do the following for each stimulus
    
    if noisefreq(1,2) == 0
        nsfreq = noisefreq(1,1) * 25;                                       % convert from Hz
    else
        nsfreq = (noisefreq(1,1) + ((rand*2*noisefreq(1,2))-noisefreq(1,2))) * 25; % randomize noise frequency
    end
    
    if noiseamp(1,2) == 0
        nsamp = noiseamp(1,1);
    else
        nsamp = noiseamp(1,1) + ((rand*2*noiseamp(1,2))-noiseamp(1,2));     % randomize noise amp
    end
    
    nsphase = rand;                                                         % randomize noise phase
    
    bslblock_spikes = mfr * (norm(bslblock)/1000);
    
    bsldata = MCsim(@(xdata) max(0,1+sin(((xdata+(nsphase*10000))*((nsfreq*pi)/10000)))*(nsamp)),1,bslblock,1,[1,bslblock],bslblock_spikes,0,0,max(1,abs(amp)+1),ylims);
    
    if a <= respstims       % if this is a trial during which the neuron is showing a response
        
        rspdata{a,1} = cell2mat(MCsim(@(xdata) max(0,(((exp(-(xdata - rspdealy).^2/(2*sigma^2)))*abs(amp))+1) + sin(((xdata+(nsphase*10000))*((nsfreq*pi)/10000)))*(nsamp)),1,10000,1,bslrange,bsl_spikes,doplot,0,max(1,abs(amp)+1),ylims));
        grndtruth(a,:) = max(0,(((exp(-(xdata - rspdealy).^2/(2*sigma^2)))*abs(amp))+1) + sin(((xdata+(nsphase*10000))*((nsfreq*pi)/10000)))*(nsamp));
        
    elseif percent < 0      % if this is a trial during which the neuron is showing a negative response
        
        namp = -(abs(amp));
        rspdata{a,1} = cell2mat(MCsim(@(xdata) max(0,(((exp(-(xdata - rspdealy).^2/(2*sigma^2)))*namp)+1) + sin(((xdata+(nsphase*10000))*((nsfreq*pi)/10000)))*(nsamp)),1,10000,1,bslrange,bsl_spikes,doplot,0,max(1,namp+1),ylims));
        grndtruth(a,:) = max(0,(((exp(-(xdata - rspdealy).^2/(2*sigma^2)))*namp)+1) + sin(((xdata+(nsphase*10000))*((nsfreq*pi)/10000)))*(nsamp));
        
    else                    % if this is a trial during which the neuron is not showing a response
        
        rspdata{a,1} = cell2mat(MCsim(@(xdata) max(0,1+sin(((xdata+(nsphase*10000))*((nsfreq*pi)/10000)))*(nsamp)),1,10000,1,bslrange,bsl_spikes,doplot,0,max(1,abs(amp)+1),ylims));
        grndtruth(a,:) = max(0,ones(1,10000)+sin(((xdata+(nsphase*10000))*((nsfreq*pi)/10000)))*(nsamp));
        
    end
    
    % Concatenate your trials into the relevant variables
    rawdata = sort([rawdata,bsldata{1,1} + (block*a),(rspdata{a,1}+bslblock) + (block*a)]);   
    peristim{a,1} = sort(rspdata{a,1} - (xdata(1,end)/2));
    peristim_cntr{a,1} = sort(((cell2mat(MCsim(@(xdata) max(0,1+sin(((xdata+(nsphase*10000))*((nsfreq*pi)/10000)))*(nsamp)),1,10000,1,bslrange,bsl_spikes,doplot,0,max(1,abs(amp)+1),ylims))) - (xdata(1,end)/2)));
    grndtruth_cntr(a,:) = max(0,1+sin(((xdata+(nsphase*10000))*((nsfreq*pi)/10000)))*(nsamp));
    stimtime(a,1) = (block * a) + (xdata(1,end)/2);
    
end

%% Add some baseline firing without response to the beginning and end of your whole experiment to pad it out

for b = 1:2
    
    if noisefreq(1,2) == 0
        nsfreq = noisefreq(1,1) * 25;                                       % convert from Hz
    else
        nsfreq = (noisefreq(1,1) + ((rand*2*noisefreq(1,2))-noisefreq(1,2))) * 25; % randomize noise frequency
    end
    
    if noiseamp(1,2) == 0
        nsamp = noiseamp(1,1);
    else
        nsamp = noiseamp(1,1) + ((rand*2*noiseamp(1,2))-noiseamp(1,2));     % randomize noise amp
    end
    
    nsphase = rand;                                                         % randomize noise phase
    
    switch b
        case 1
            bslblock_spikes = mfr/(1000/norm(block));
            bsldata = MCsim(@(xdata) max(0,1+sin(((xdata+(nsphase*10000))*((nsfreq*pi)/10000)))*(nsamp)),1,block,1,[1,block],bslblock_spikes,0,0,max(1,abs(amp)+1),ylims);
        case 2
            rawdata = [rawdata,bsldata{1,1} + (block*(a+1))];
            bsldata = MCsim(@(xdata) max(0,1+sin(((xdata+(nsphase*10000))*((nsfreq*pi)/10000)))*(nsamp)),1,block,1,[1,block],bslblock_spikes,0,0,max(1,abs(amp)+1),ylims);
            rawdata = [bsldata{1,1},rawdata];
    end
    
end

end


function [data] = MCsim(f,xMin,xMax,ntrials,bslrange,bsl_spikes,plotit,yMin,yMax,ylims)

close all
data{ntrials,1} = [];
refract = 3;    % Refractory period for your neuron in ms (set to zero if you don't wan't a refractory period)

for a = 1:ntrials
    
    xValues = 0;
    yValues = 0;
    
    if plotit   % plot your MC Simulation
        figure('color', 'white');
        title('Monte Carlo simulation in progress','FontSize',18)
        ylim(ylims);
        hold on
        delta = (xMax - xMin) / 1e4;
        x = xMin: delta : xMax;
        y = f(x);
        plot(x,y,'r','LineWidth',2)
        curve = plot(xValues, yValues, 'ok','MarkerSize',7,'MarkerFace',[0 0 0]); % 'lineWidth', 1.25);
        set(curve, 'XDataSource', 'xValues');
        set(curve, 'YDataSource', 'yValues');
        xlabel('You can disable these plots by setting the variable "plotit" to zero','FontSize',16)
    end
    
    %% Generate you MC data
    
    getstarted = 1;
    
    while getstarted == 1 || size(xValues((xValues(1,:) > bslrange(1,1)) & (xValues(1,:) < bslrange(1,2))),2) < bsl_spikes  % generate spikes until needed number is reached
        
        x = xMax - (xMax - xMin)*rand();        % generate x between xMin and xMax
        y = yMax - (yMax - yMin)*rand();        % generate y between 0 and yMax
        
        if y <= f(x) && (x - max(xValues(1,xValues(1,:) < x))) > refract   % Accept a spiek only if y lies below (or on) the curve & isi is > refract
            
            xValues(1,end+1) = x;
            yValues(1,end+1) = y;
            
            if plotit
                refreshdata(curve, 'caller')        % update the plot
                drawnow
            end
            
        end
        getstarted = 0;
    end
    
    data{a,:} = xValues;
    
end

close all

end



%% Remnants
% 
% %% Recreate the cluster_class file format
% 
%  cluster_class(:,2) = rawdata;
%  cluster_class(:,1) = 1;
% 
% %% Recreate the rawstims file format
% 
% rawstims(:,2) = stimtime;
% rawstims(:,1) = 1;
% rawstims = sortrows(rawstims,2);

