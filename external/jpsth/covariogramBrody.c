/*
	mex code for Covariogram in Brody 1999
*/

/*
	covariogramBrody.c

	Copyright 2008 Vanderbilt University.  All rights reserved.
	John Haitas, Jeremiah Cohen, and Jeff Schall

	This file is part of JPSTH Toolbox.

	JPSTH Toolbox is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	JPSTH Toolbox is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with JPSTH Toolbox.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "math.h"
#include "mex.h"
#include "matrix.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int i, j, k, start, finish, lag, currentLag, trialLength, trials;
	double *n_1, *n_2, *p_1, *p_2, *s_1, *s_2;
	double crossCorr, shuffleCorrector, s1s2, p1s2, s1p2, sigma;
	double *covariogram, *sigHigh, *sigLow;
	
	/* get dimensions of spike counts */
	trials = mxGetM(prhs[0]);
	trialLength = mxGetN(prhs[0]);
		
	/* connect input arguments for MATLAB interface */
	n_1 = mxGetData(prhs[0]);
	n_2 = mxGetData(prhs[1]); 
	p_1 = mxGetData(prhs[2]);
	p_2 = mxGetData(prhs[3]);
	s_1 = mxGetData(prhs[4]);
	s_2 = mxGetData(prhs[5]);
	
	if (nrhs != 7) { lag = 50; }
	else { lag = mxGetScalar(prhs[6]); }
	
	/* connect output to MATLAB interface */
	plhs[0] = mxCreateDoubleMatrix(1,2*lag+1,mxREAL);
	covariogram = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(1,2*lag+1,mxREAL);
	sigHigh = mxGetPr(plhs[1]);
	plhs[2] = mxCreateDoubleMatrix(1,2*lag+1,mxREAL);
	sigLow = mxGetPr(plhs[2]);
		
	/* algorithm for Covariogram in Brody 1999 */
	for (i = 0; i < 2*lag+1; i++) {
		currentLag = i - lag;
		start = (currentLag < 0) ? -currentLag : 0;
		finish = (currentLag < 0) ? trialLength : trialLength - currentLag;
		crossCorr = 0.;
		shuffleCorrector = 0.;
		s1s2 = 0;
		p1s2 = 0;
		s1p2 = 0;
		for (j = start; j < finish; j++) {
			for (k = 0; k < trials; k++) 
				crossCorr += n_1[k + j*trials] * n_2[k + (j+currentLag)*trials];
			shuffleCorrector += ( p_1[j] * p_2[j+currentLag] );
			s1s2 += ( pow(s_1[j],2) * pow(s_2[j+currentLag],2) );
			p1s2 += ( pow(p_1[j],2) * pow(s_2[j+currentLag],2) );
			s1p2 += ( pow(s_1[j],2) * pow(p_2[j+currentLag],2) );
		}
		sigma = sqrt((s1s2+p1s2+s1p2)/trials);
		
		covariogram[i] = (crossCorr / trials) - shuffleCorrector;
		sigHigh[i] = 2 * sigma;
		sigLow[i] = -2 * sigma;
	}

}
