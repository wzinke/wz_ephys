% COVARIOGRAMBRODY computes a covariogram of spike_1 and spike_2 at a 
% given lag
%
% to compile enter the following command at the MATLAB Command Window
% >> mex('covariogramBrody.c')
%
% MATLAB equivalent to covariogramBrody.c listed below
%

% covariogramBrody.m
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

function [thisCovariogram, sigHigh, sigLow] = covariogramBrody(spike_1, spike_2, p1, p2, s1, s2, lag)
    % warn user that they should use compiled version
    warning('JPSTH:compileMex', ...
    		['You are using the MATLAB version of covariogramBrody()\n' ...
				'For better performance please compile covariogramBrody.c with this command: mex(\''covariogramBrody.c\'')']);

	if nargin < 7 lag = 50; end
	
	trials = size(spike_1,1);
	trialLength =  size(spike_1,2);
		
	s1s2 = zeros(1,2*lag+1);
	p1s2 = zeros(1,2*lag+1);
	s1p2 = zeros(1,2*lag+1);
	crossCorr =  zeros(1, 2*lag+1);
	shuffleCorrector = zeros(1,2*lag+1);
	
	for i = 1:2*lag+1
		currentLag = i - lag - 1;
		
		if currentLag < 0 
			jVector = 1-currentLag:trialLength;
		else 
			jVector = 1:trialLength-currentLag; 
		end
		
		for j = jVector
			crossCorr(i) = crossCorr(i) + mean(spike_1(:,j) .* spike_2(:,j+currentLag));
			shuffleCorrector(i) = shuffleCorrector(i) + p1(j) * p2(j+currentLag);
			s1s2(i) = s1s2(i) + s1(j).^2 * s2(j+currentLag).^2;
			p1s2(i) = p1s2(i) + p1(j).^2 * s2(j+currentLag).^2;
			s1p2(i) = s1p2(i) + s1(j).^2 * p2(j+currentLag).^2;
		end
	end
	
	thisCovariogram = crossCorr - shuffleCorrector;
	
	sigma = sqrt((s1s2 + p1s2 + s1p2)/trials);
	sigHigh = 2 * sigma;
	sigLow = -2 * sigma;
end
