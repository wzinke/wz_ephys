% SPIKE_COUNTS converts timeStamps into logical time-indexed spikes
% returns makes matrix of ones and zeros
% rows are trials
% columns are times

% spikeCounts.m
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

function [thisSpikeCount] = spikeCounts(timeStamps, timeWindow, binWidth)

	if nargin < 2, timeWindow = [-100 500]; end
	if nargin < 3, binWidth = 1; end
		
	binCenters = timeWindow(1)+(binWidth/2):binWidth:timeWindow(2)-(binWidth/2);
	
	thisSpikeCount = hist(timeStamps',binCenters)';
end
