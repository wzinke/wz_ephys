function [index] = cumulativeSaccadePotential(dise_time, min_dise, max_dise)
%Function [index] = cumulativeSaccadePotential(dise_time, min_dise, max_dise)
%
% Function calculates the cumulative saccade potential index, see Leppanen,
% et al. 2014.

index = 1 - ((max_dise - dise_time)/(max_dise - min_dise));