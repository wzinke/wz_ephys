function ff = SPK_get_Fano(spikecounts)
% wz_spk_get_Fano - calculate fano factor for spike times
%  
% DESCRIPTION 
% 	This function calculates the fano factor, i.e. the variation of spike
% 	counts across trials divided by the mean of spike counts across trials.
%  
% SYNTAX 
% ff = wz_spk_get_Fano(spikecounts)

%   Input:
%       spikecounts  vector of spike counts or matrix with spike times.
%                    Missing values must be NaN!
%
%   Output:
%       ff    fano factor
%
% REFERENCES 
%
% [-] Fano U (1947)
%     Ionization yield of radiations. II. The fluctuations of the number of ions.
%     Phys Rev 72(1), 26?9.
%
% [-] Gur M, Beylin A, Snodderly DM (1997)
%    Response variability of neurons in primary visual cortex (V1) of alert monkeys.
%    J Neurosci 17(8), 2914?20.
%
% ......................................................................... 
% wolf zinke, wolfzinke@gmail.com 
%
% $Created : 30-Jul-2014 by wolf zinke
% $Modified: 

% ========================================================================= 
%% check input arguments and set default values 
if(min(size(spikecounts)) > 1)  % input is a 2d matrix with spike times
    spikecounts = sum(isfinite(spikecounts),2);
end

% ========================================================================= 
%% prepare input data &  do the calculation

if(length(spikecounts) > 4)
    ff = var(spikecounts) / mean(spikecounts);
else
    ff = NaN;
end

