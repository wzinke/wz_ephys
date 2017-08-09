function lv = SPK_get_ISIlv(isi, do_trial)
% wz_spk_get_ISIcv - coefficient of variation for an ISI distribution
%  
% DESCRIPTION 
% Calculation of the local variation of ISIs as introduced by Shinomoto et al.
% Here a value is calculated for two consecutive ISIs, where LV is the mean of 
% all these values. Even though it wasn not decribed by the authors this way,
% this function also has the time resolved LV-values as additional output.
%
% Times are specified as shown here:
% S1-3: occurrence of spike
%         S1     S2        S3
%       __|______|_________|__
%         |ISI(i)| ISI(i+1)|
%            Time of tIV
%
% SYNTAX 
% cv = wz_spk_get_ISIlv(isi, do_trial)

%   Input:
%       isi       vector or matrix containin ISI's. Missing values must be NaN!
%
%       do_trial  flag that specifies if calculation should be done per row
%
%   Output:
%       lv    index of local variation
%
% REFERENCES 
%
%   [1] Shinomoto S, Shima K, and Tanji J (2003).
%       Differences in Spiking Patterns among Cortical Neurons.
%       Neural Computation, 15, 2823-2842
%
%   [2] Shinomoto S, Miura K, Koyama S (2005).
%       A Measure of Local Variation of Inter-Spike Intervals.
%       Biosystems,79(1-3), 67-72
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

% ========================================================================= 
%% do the calculation
if(sum(sum(isfinite(isi),2)) > 4)
    
    % calculate lv for single pairs
    lvsgl =    3 * (isi(:,2:end) - isi(:,1:end-1)).^2 ./ ...
                   (isi(:,2:end) + isi(:,1:end-1)).^2 ;
                     
    % get the average cv2 value
    if(do_trial == 0)
        lvsgl = lvsgl(:);
        lvsgl(~isfinite(lvsgl)) = [];
        
        lv =    mean(lvsgl);
    else
        lv = nanmean(lvsgl, 2);
    end
else  % not sufficient data to get a meaningful value
    if(do_trial == 0)
        lv = NaN;
    else
        lv = nan(size(isi,1),1);
    end
end

