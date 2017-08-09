function [spiketimes, ordr] = wz_spk_clip(spiketimes,clip)
% clip spikes that occur after a specified event.
%  
% DESCRIPTION 
% This function removes spike times that occur after an event at the time
% specified with <clip>. This clipping could be done on a trial basis.
% 
% SYNTAX 
% SPKobj = wz_get_SPKobj(spk_times, tmwin, tzero, clip, na_val, prec, reso)
%
%   Input:
%         <spk_times> should be a 2D matrix containing spike timies with rows representing
% ......................................................................... 
% wolf zinke, wolfzinke@gmail.com 

nTrials = size(spiketimes,1);

if(length(clip) == 1)
    spiketimes(spiketimes>clip) =  NaN;
    ordr = 1:nTrials;
elseif(length(clip) == nTrials)
    [~ , ordr] = sort(clip,'descend');
end
spiketimes(bsxfun(@ge, spiketimes, clip(:))) = NaN;

