% SIGNIFICANT_SPAN identifies longest timespan that exceeds significance
% if one exists

% significantSpan.m
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

function [spanEndpoints] = significantSpan(vectorA, sig)
	if length(sig)==1
		sig = repmat(sig,1,length(vectorA));
	end
		
	% find indices that exceed significance
	difference = vectorA > sig;

    % after this you will see -1 on the index before a span starts
    % you will see a 1 on the index where a span ends
    divDifference = [difference(2:end) 0] - difference(:)';
        
    % find all the beginning indices
    beginnings = find(divDifference==1) + 1;
        
    % find all the ending indices
    endings = find(divDifference==-1);

	% if the first index is significant beginnings wouldn't have caught it  	 	 
    if ~isempty(endings) && length(beginnings) < length(endings) 		 
        beginnings = [1 beginnings]; 		 
    end

    % identify the longest span(s)
	maxIndices = (endings-beginnings)==max(endings-beginnings);
	
    %find the indices of the longest spans
    spanEndpoints = [];
    theIndices = find(maxIndices==1);
    for i=1:length(theIndices)
        theIndex = theIndices(i);
        spanEndpoints = horzcat(spanEndpoints,[beginnings(theIndex) endings(theIndex)]);
    end
end
