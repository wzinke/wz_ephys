function cv = SPK_get_ISIcv(isi, do_trial)
% wz_spk_get_ISIcv - coefficient of variation for an ISI distribution
%  
% DESCRIPTION 
% 	This function calculates the coefficient for the interspike interval
% 	distribution
%  
% SYNTAX 
% cv = wz_spk_get_ISIcv(isi, do_trial)

%   Input:
%       isi       vector or matrix containin ISI's. Missing values must be NaN!
%
%       do_trial  flag that specifies if calculation should be done per row
%
%   Output:
%       cv    coeeficient of variation
%
% REFERENCES 
%
% [-] Softky WR & Koch C (1993)
%     The High Irregular Firing of Cortical Cells is Inconsistent with 
%     Temporal Integration of Random EPSPs.
%     J Neurosci 13(1), 334-350
%
% ......................................................................... 
% wolf zinke, wolfzinke@gmail.com 
%
% $Created : 30-Jul-2014 by wolf zinke
% $Modified: 

% ========================================================================= 
%% check input arguments and set default values 
if(~exist('do_trial','var')   || isempty(do_trial))
    do_trial = 0;
end

% ========================================================================= 
%% prepare input data &  do the calculation

if(do_trial == 0)
    isi = isi(:);          % make sure the data is a single vector
    isi(isnan(isi)) = [];  % remove NaN's from the data
    cv  = std(isi) ./ mean(isi);
else
    cv  = nanstd(isi,[],2) ./ nanmean(isi,2);
end


