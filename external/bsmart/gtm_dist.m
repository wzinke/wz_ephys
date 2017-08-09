function [DIST, minDist, maxDist] = gtm_dist(T, Y, m)
% Calculate the squared distances between two sets of data points. 
%
%		This function calculates distances between all data 
%		points in the two data sets T and Y and returns
%		them in a matrix.
%
% Synopsis:	[DIST, minDist, maxDist] = gtm_dist(T, Y, m)
% 		[DIST] = gtm_dist(T, Y)
%
% Arguments:	T, Y -	data set matrices in which each row is a
%			data point; dimensions N-by-D and K-by-D
%			respectively
%
%		m - 	mode of calculation; iff m > 0, min- and  
%			maxDist (below) are calculated; the 
%			default mode is 0 
%
% Return:	DIST -	matrix containing the calculated distances;
%			dimension K-by-N; DIST(k,n) contains the
%			squared distance between T(n,:) and Y(k,:).
%
%		minDist,
%		  maxDist -	vectors containing the minimum and 
%				maximum of each column in DIST, 
%				respectively; 1-by-N; required 
%				iff m > 0.
%			
% Notes:	This m-file provides this help comment and a MATLAB
%		implementation of the distance calculation. If, 
%		however, a mex-file with the same name is present
%		in the MATLABPATH, this will be called for doing
%		the calculation. As this is a computationally 
%		demanding step of the algorithm, an efficient
%		mex-file implementation will improve the performance 
%		of the GTM training algorithm.
%
% See also:	gtm_dstg
%

% Version:	The GTM Toolbox v1.0 beta
%
% Copyright:	The GTM Toolbox is distributed under the GNU General Public 
%		Licence (version 2 or later); please refer to the file 
%		licence.txt, included with the GTM Toolbox, for details.
%
%		(C) Copyright Markus Svensen, 1996

if (nargin == 3)
  if (m < 0)
    error(['Invalid value for mode: ', num2str(m)]);
  else
    mode = m;
  end
elseif (nargin == 2)
  mode = 0;
else
  error('Wrong number of input arguments');
end

if (mode > 0)
  if (nargout < 3)
    error('Calculation mode > 0 requires 3 output arguments!');
  end
end

[N, tD] = size(T);
[K, yD] = size(Y);

if (yD ~= tD)
  error('Mismatch in number of columns between T and Y.');
end

% Summing over components of matrices, we can make use of a 'matrix 
% version' of the identity: (a - b)^2 == a^2 + b^2 - 2*a*b, 
% which yields faster execution. When T and Y consist of single columns 
% (which may be the case with calls from gtm_gbf), this has to be handled 
% as a special case. 
if yD > 1
  DIST = Y*T';
  tt = sum(T'.^2);
  yy = sum(Y'.^2);
  DIST = yy'*ones(1,N) + ones(K,1)*tt - 2*DIST;
else
  DIST = zeros(K, N);
  for n=1:N
    for k=1:K
      DIST(k,n) = sum((T(n,:) - Y(k,:)).^2);
    end
  end
end

if (mode > 0) 
  minDist = min(DIST);
  maxDist = max(DIST);
end
