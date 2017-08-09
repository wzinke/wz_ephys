/*	equation3.c

	mex code for Eq.3 in Aertsen et al. 1989
*/

/*
	equation3.c
	
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
	int u, v, k, windowSize, trials;
	double *n_1, *n_2, *rawJPSTH, sum;
	
	/* get dimensions of spike counts */
	trials = mxGetM(prhs[0]);
	windowSize = mxGetN(prhs[0]);
	
	/* connect input arguments for MATLAB interface */
	n_1 = mxGetData(prhs[0]);
	n_2 = mxGetData(prhs[1]);
	
	/* connect output to MATLAB interface */
	plhs[0] = mxCreateDoubleMatrix(windowSize,windowSize,mxREAL);
	rawJPSTH = mxGetPr(plhs[0]);
	
	/* Algorithm for Equation 3 in Aertsen et al. 1989 */
	for (u = 0; u < windowSize; u++) {
		for (v = 0; v < windowSize; v++) {
			sum = 0.;
			for (k = 0; k < trials; k++) 
				sum += n_1[k + u*trials] * n_2[k + v*trials];
			rawJPSTH[u + v*windowSize] = sum / trials;
		}
	}
}
