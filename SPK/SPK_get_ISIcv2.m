function cv2 = SPK_get_ISIcv2(isi, do_trial)
% wz_spk_get_ISIcv - coefficient of variation for adjacent ISI pairs
%  
% DESCRIPTION 
% To measure the intrinsic variability of spiketrains Holt et al. introduced a
% new  measure (CV2). This value is more robust agains changes of the firingrate, 
% because it uses two adjacent ISIs for the calculation. In words, CV2is calculated 
% by deviding the standard deviation of two consecutive ISIsby their mean and 
% multiplying this by with sqrt(2) to get cv2 = 1 for a Poisson process:
%
%           CV2 = (2|ISI(i+1) - ISI(i)|) / (ISI(i+1) + ISI(i))
%
% Times are specified as shown here:
% S1-3: occurrence of spike
%         S1      S2        S3
%       __|_______|_________|__
%         |ISI(i) | ISI(i+1)|
%            Time of tIV
%  
% SYNTAX 
% cv = wz_spk_get_ISIcv2(isi, do_trial)

%   Input:
%       isi       vector or matrix containin ISI's. Missing values must be NaN!
%
%       do_trial  flag that specifies if calculation should be done per row
%
%   Output:
%       cv2    coeeficient of variation for adjacent ISI pairs
%
% REFERENCES 
%
% [-] Softky WR & Koch C (1993)
%     The High Irregular Firing of Cortical Cells is Inconsistent with 
%     Temporal Integration of Random EPSPs.
%     J Neurosci 13(1), 334-350
%
% [-] Holt GR, Softky WR, Koch C, and Douglas RJ (1996).
%     Comparison of Discharge Variability In Vitro and In Vivo in Cat Visual Neurons.
%     J Neurophysiol, 75(5), 1806-1814
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
    
    % calculate cv2 for single pairs
    cv2sgl = 2 * abs(isi(:,2:end) - isi(:,1:end-1)) ./ ...
                    (isi(:,2:end) + isi(:,1:end-1));
                 
    % get the average cv2 value
    if(do_trial == 0)
        cv2sgl = cv2sgl(:);
        cv2sgl(~isfinite(cv2sgl)) = [];
        
        cv2 =    mean(cv2sgl);
    else
        cv2 = nanmean(cv2sgl, 2);
    end
else  % not sufficient data to get a meaningful value
    if(do_trial == 0)
        cv2 = NaN;
    else
        cv2 = nan(size(isi,1),1);
    end
end

