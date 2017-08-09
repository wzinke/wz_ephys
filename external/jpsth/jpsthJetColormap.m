% JPSTHJETCOLORMAP produce a colormap with warm colors for positive
% values and cold colors for negative values with green representing
% zero

% jpsthJetColormap.m
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

function [thisColormap] = jpsthJetColormap(cmLength)
	if nargin < 1, cmLength=65; end
	if mod(cmLength,2)==0, cmLength=cmLength+1; end
	thisColormap = zeros(cmLength, 3);
	
	clim = get(gca,'CLim');
	cMin = clim(1);
	cMax = clim(2);
	zeroI = fix((0-cMin)/(cMax-cMin)*cmLength)+1;
	
	negLength = zeroI - 1;
	posLength = cmLength - zeroI;
	
	for i = 1:2*negLength/3
		j = i-1;
		thisColormap(i,:) = [0 j/(2*negLength/3) 1];
	end
	
	startIndex = fix(2*negLength/3)+1;
	for i = startIndex:zeroI-1
		j = length(startIndex:zeroI-1) - (i - startIndex) - 1;
		thisColormap(i,:) = [0 1 (3*j/negLength)^.7];
	end
	
	thisColormap(zeroI,:) = [0 1 0];
		 
	startIndex = zeroI+1;
	for i = startIndex:zeroI+(posLength/3)
		j = i - startIndex;
		thisColormap(i,:) = [(3*j/(posLength))^.7 1 0];
	end
	
	startIndex = fix(zeroI+(posLength/3))+1;
	for i = startIndex:cmLength
		j = length(startIndex:cmLength) - (i - startIndex) - 1;
		thisColormap(i,:) = [1 (j/(2*posLength/3)) 0];
	end
	
end
