% ALIGNTIMESTAMPS align signal by subtracting event timstamp from spike
% timeStamps

% alignTimeStamps.m
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

function [ alignedTimeStamps ] = alignTimeStamps( timeStamps, event)
	% set 0's and spikes outside window to NaN
	timeStamps(timeStamps==0) = NaN;
	
	% subtract event time to make all timeStamps relative to this event
	alignedTimeStamps = timeStamps-repmat(event(:),1,size(timeStamps,2));
	
end
