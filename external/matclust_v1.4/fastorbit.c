

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

	mxLogical 			  *logicptr;
	int               i, j, count, memvector, remainder, tmpbit;
        int            tempnum;
        int           *input;
	double            *bits;
	mxLogical	          *output;
	int				  ninput;
        int                               inputlength, ncolumns;
	int				  nbits;


    /* Check numbers of arguments */
    if (!(nrhs == 2) || !(nlhs == 1)){
        Usage();
    }
   





	input = mxGetData(prhs[0]);
	ninput = mxGetM(prhs[0]) * mxGetN(prhs[0]);
        inputlength = mxGetM(prhs[0]);
        ncolumns = mxGetN(prhs[0]);

	bits = mxGetPr(prhs[1]);
	nbits = mxGetM(prhs[1]) * mxGetN(prhs[1]);
        

	
	output = (mxLogical *) mxCalloc(inputlength, sizeof(mxLogical));

     
         if (!mxIsInt32(prhs[0])) {
           Usage();
        }           
        
        
	for (i = 0; i < inputlength; i++) {
		tempnum = 0;
                for (j = 0; j < nbits; j++) {
                    if (!tempnum) {
                       tmpbit = (int) bits[j] - 1;
                       memvector = tmpbit>>5;    /*equiv to floor(tmpbit/32)*/                 
                       remainder = (tmpbit) & 31; /* find the modulus of 32*/ 
                      
                        
                       tempnum = tempnum || (( (input[memvector*inputlength + i])>>remainder) & 01);
                       /*mexPrintf("vector %d remainder %d %d\n",memvector, remainder, tempnum);*/
                    }
                    else {
                       break;
                    }
                }
                
                output[i] = (mxLogical) tempnum;
        }

	plhs[0] = mxCreateLogicalMatrix(inputlength, 1);

	logicptr = mxGetData(plhs[0]);
	/* copy into the output array */
	for (i = 0; i < inputlength; i++) {
		 *logicptr = output[i];
		 logicptr++;
	}


	mxFree(output);

	return;
}

void Usage(void)
{
	mexErrMsgTxt("wrong usage\n");
}


