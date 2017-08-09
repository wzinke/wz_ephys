function ir = SPK_get_ISIir(isi, do_trial)
% wz_spk_get_ISIcv - coefficient of variation for an ISI distribution
%  
% DESCRIPTION 
% based on a new method introduced by SN Baker this function serves to measure
% time dependent changes in the irregularity of neural spiking.
%
% This measure is based on a simple calculation IR = abs(log(ISI(i) / ISI(i+1))),
% where ISI(i) and ISI(i+1) are two subsequent inter spike intervals.
% The corresponding time for each value should be the time of the centre
% spike defining this two intervals.

% Times are specified as shown here:
% S1-3: occurrence of spike
%         S1      S2        S3
%       __|______|_________|__
%         |ISI(i)| ISI(i+1)|
%            Time of tIV

% To correct for artefacts induced by fast rate changes it is possible to use
% the signed values sIR = log(ISIa / ISIb). Usually the histogram of the sIR
% values should be a gamma distributuion centered at 0. When a artefact occurs
% it is possible, to remove values beyond a specified threshold and then to
% calculate the IR with the remaining data.
%  
% SYNTAX 
% cv = wz_spk_get_ISIir(isi, do_trial)

%   Input:
%       isi       vector or matrix containin ISI's. Missing values must be NaN!
%
%       do_trial  flag that specifies if calculation should be done per row
%
%   Output:
%       ir    spiking irregularity index
%
% REFERENCES 
%
% [1] Davies RM, Gerstein GL, and Baker SN (2006).
%     Measurement of Time-Dependent Changes in the Irregularity of Neural Spiking.
%     J Neurophysiol, 96(2), 906-18.
%
%
% ......................................................................... 
% wolf zinke, wolfzinke@gmail.com 
%
% $Created : 30-Jul-2014 by wolf zinke
% $Modified: 

% ========================================================================= 
%% check input arguments and set default values 
if(~exist('do_trial','var') || isempty(do_trial))
    do_trial = 0;
end

% correct for ISIs of duration 0 ms
zcnt = sum(isi(:) == 0);

if(zcnt > 0)
    warning([int2str(zcnt), ' ISIs with 0 ms duration found and removed! [',num2str(zcnt/sum(isfinite(isi(:))),'%.4f'),'%]']);
    isi = wz_cropNaN(isi, [], [], 0);
end

% ========================================================================= 
%% do the calculation
if(max(sum(isfinite(isi),2)) > 4)
    
    % calculate ir for single pairs
    irsgl =   abs(log(isi(:,2:end) ./ isi(:,1:end-1)));
                     
    % get the average ir value
    if(do_trial == 0)
        irsgl = irsgl(:);
        irsgl(~isfinite(irsgl)) = [];
        
        ir =    mean(irsgl);
    else
        ir = nanmean(irsgl, 2);
    end
else  % not sufficient data to get a meaningful value
    if(do_trial == 0)
        ir = NaN;
    else
        ir = nan(size(isi,1),1);
    end
end

