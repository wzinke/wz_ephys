% EQUATION3 computes the raw JPSTH of spike_1 and spike_2
%
% to compile enter the following command at the MATLAB Command Window
% >> mex('equation3.c')
%
% MATLAB equivalent to equation3.c listed below

% equation3.m
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

function [rawJPSTH] = equation3(spike_1,spike_2)
    % warn user that they should use compiled version
    warning('JPSTH:compileMex', ...
    		['You are using the MATLAB version of equation3()\n' ...
				'For better performance please compile equation3.c with this command: mex(\''equation3.c\'')']);
    
	rawJPSTH = zeros(size(spike_1,2));							
	for u = 1:size(spike_1,2)										
		for v = 1:size(spike_2,2)									
			rawJPSTH(u,v) = mean(spike_1(:,u) .* spike_2(:,v));		
		end														
	end															
end
