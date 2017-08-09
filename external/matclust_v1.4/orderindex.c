

#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>



void Usage(void) ;


/******************************************************************************
  INTERFACE FUNCTION
  */

void mexFunction(
    int           nlhs,           /* number of expected outputs */
    mxArray       *plhs[],        /* array of pointers to output arguments */
    int           nrhs,           /* number of inputs */
    const mxArray *prhs[])         /* array of pointers to input arguments */
{

	double 	      *dblptr;
	int               i, j, count;
        int            tempnum;
        double           *index;
	double            *clusterson;
	double	          *output;
	int				  nindex, indexlength, indexwidth, nclusters;


    /* Check numbers of arguments */
    if (!(nrhs == 2) || !(nlhs == 1)){
        Usage();
    }
   





	index = mxGetPr(prhs[0]);
	nindex = mxGetM(prhs[0]) * mxGetN(prhs[0]);
        indexlength = mxGetM(prhs[0]);
        indexwidth = mxGetN(prhs[0]);

	clusterson = mxGetPr(prhs[1]);
	nclusters = mxGetM(prhs[1]) * mxGetN(prhs[1]);
        

	
	output = (double *) mxCalloc(nindex, sizeof(double));

     
       count = 0;
        
	for (i = 1; i < nindex; i+=2) {
		if (index[i] == 2) {
                   output[count] = index[i-1];
                   output[count+1] = index[i];
                   count+=2;
                   
                }
        }
        
        /*for (i = 1; i < nindex; i+=2) {
		if (index[i] == 1) {
                   output[count] = index[i-1];
                   output[count+1] = index[i];
                   count+=2;
                   
                }
        }       */
        for (j = 0;j<nclusters;j++) {
           for (i = 1; i < nindex; i+=2) {
		if (index[i] == clusterson[j]+2) {
                   output[count] = index[i-1];
                   output[count+1] = index[i];
                   count+=2;
                  
                }
           } 
        }      
       

	plhs[0] = mxCreateDoubleMatrix(indexlength, indexwidth, mxREAL);

	dblptr = mxGetPr(plhs[0]);
	/* copy into the output array */
	for (i = 0; i < nindex; i++) {
		 *dblptr = output[i];
		 dblptr++;
	}


	mxFree(output);

	return;
}

void Usage(void)
{
	mexErrMsgTxt("wrong usage\n");
}


