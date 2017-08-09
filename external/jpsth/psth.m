% PSTH computes the PSTH of a single spike train
% Equations from Aertsen et al. 1989
% psthMean:		Equation 1
% psthStdDev:		for Equation 7a
% psthVariance:		for Brody Covariogram

% psth.m
% 
% Copyright 2008 Vanderbilt University.  All rights reserved.
% John Haitas, Jeremiah Cohen, and Jeff Schall
% 
% This file is part of JPSTH Toolbox.
% 
% JPSTH Toolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% JPSTH Toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with JPSTH Toolbox.  If not, see <http://www.gnu.org/licenses/>.

function [psthMean, psthStdDev, psthVariance] = psth(spikeCounts)
	if sum(sum(spikeCounts))==0
		psthException = MException('PSTH:NoSpikes','There are no spikes');
		throw(psthException);
	end
	% take the mean, standard deviation, and variance from the spike counts
	psthMean = mean(spikeCounts);
	psthStdDev = std(spikeCounts);
	psthVariance = var(spikeCounts);
end
