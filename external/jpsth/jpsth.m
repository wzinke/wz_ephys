% JPSTH calculates jpsth for two spike counts aligned around an event

% jpsth.m
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

function [jpsthData] = jpsth(n_1, n_2, coinWidth)
	% Supress divideByZero warning for Eq. 9
	warning('off','MATLAB:divideByZero');

	if nargin < 3, coinWidth = 10; end

	% calculate PSTH for each signal, Eq. 1 and terms for Eq. 7a
	[psth_1, psth_1_sd, psth_1_var] = psth(n_1);
	[psth_2, psth_2_sd, psth_2_var] = psth(n_2);
	
	% JPSTH Equations from Aertsen et al. 1989
	rawJPSTH = equation3(n_1, n_2);						% Eq. 3 
	psthOuterProduct = psth_1(:) * psth_2(:)';			% Eq. 4
	unnormalizedJPSTH = rawJPSTH - psthOuterProduct;	% Eq. 5
	normalizer = psth_1_sd(:) * psth_2_sd(:)';			% Eq. 7a
	normalizedJPSTH = unnormalizedJPSTH ./ normalizer;	% Eq. 9
	
	% When normalizer = 0 it causes divide by zero which results in NaN
	% the following replaces these NaN's with 0's
	% this occurs for discrete time points with no spikes
	normalizedJPSTH(isnan(normalizedJPSTH)) = 0;
	
	% analyze JPSTH
	xcorrHist = crossCorrelationHistogram(normalizedJPSTH);
	pstch = pstCoincidenceHistogram(normalizedJPSTH,coinWidth);
	
	[covariogram,sigHigh,sigLow] = covariogramBrody(n_1,n_2,psth_1,psth_2,psth_1_var,psth_2_var);
	sigPeakEndpoints = significantSpan(covariogram,sigHigh);
	sigTroughEndpoints = significantSpan(-covariogram,-sigLow);
	
	
	% build data structure to return from function
	jpsthData = struct('psth_1',psth_1,'psth_2',psth_2, ...
						'normalizedJPSTH',normalizedJPSTH, ...
						'xcorrHist',xcorrHist, ...
						'pstch',pstch, ...
						'covariogram',covariogram, ...
						'sigLow',sigLow, ...
						'sigHigh',sigHigh, ...
						'sigPeakEndpoints',sigPeakEndpoints,...
						'sigTroughEndpoints',sigTroughEndpoints);
end	
			