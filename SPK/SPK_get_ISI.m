function [isi, isimat] = SPK_get_ISI(spktimes)
% wz_spk_get_ISI - calculate interspike intervals from spike times
%  
% DESCRIPTION 
% 	This function calculates interspike intervals for 
% 	
%  
% SYNTAX 
% cv = wz_spk_get_ISIcv(isi)

%   Input:
%       spktm   2D matrix containing spiketimes, with rows corresponding to trials.
%               Missing values must be NaN!
%
%   Output:
%       cv     coeeficient of variation
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
% ToDo:  - implement order of ISI to get values not to next spike, but the
%          spike following after one, two , or more intermediate spikes. 
%
% $Created : 30-Jul-2014 by wolf zinke
% $Modified: 

%  ======================================================================== 
%% prepare input data 

spktimes  = wz_cropNaN(spktimes);  % remove columns containing only NaN's

%  ======================================================================== 
%% do the calculation

isimat = diff(spktimes,1,2);

isi = isimat(:);
isi(isnan(isi)) = []; % remove NaN's from the data



