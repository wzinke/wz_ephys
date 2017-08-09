

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

	int 			  *intptr;
	int               i, j, count, memvector, remainder, tmpbit, addbit, tmpbitchange;
        int            tempnum;
        int           *input;
	double            *bits;
        mxLogical            *bitval;
	int	          *output;
	int				  ninput, nbitval;

	int				  nbits;


    /* Check numbers of arguments */
    if (!(nrhs == 3) || !(nlhs == 1)){
        Usage();
    }
   





	input = mxGetData(prhs[0]);
	ninput = mxGetM(prhs[0]) * mxGetN(prhs[0]);
        
	bits = mxGetPr(prhs[1]);
	nbits = mxGetM(prhs[1]) * mxGetN(prhs[1]);
        
        bitval = mxGetLogicals(prhs[2]);
	nbitval = mxGetM(prhs[2]) * mxGetN(prhs[2]);
        

	
	output = (int *) mxCalloc(ninput, sizeof(int));

             
        addbit = (0 | 01);
        addbit = addbit<<((int) bits[0]);
        
        if (((nbitval!=1)&&(nbitval!=ninput)) || (!mxIsLogical(prhs[2])) ){
           Usage();
        }
	for (i = 0; i < ninput; i++) {
		
                if (nbitval > 1) {
                   
                   tmpbit =  (((input[i] >> (((int) bits[0])-1)) & 01) ^ (((int) bitval[i]) & 01)) << (((int) bits[0])-1);
                   
                   
                   output[i] = (input[i] ^ tmpbit);
                }
                else {
                   tmpbit =  (((input[i] >> (((int) bits[0])-1)) & 01) ^ (((int) bitval[0]) & 01)) << (((int) bits[0])-1);
                   
                   output[i] = (input[i] ^ tmpbit);
                } 
        }

	plhs[0] = mxCreateNumericMatrix(ninput, 1, mxINT32_CLASS, mxREAL);

	intptr = mxGetData(plhs[0]);
	/* copy into the output array */
	for (i = 0; i < ninput; i++) {
		 *intptr = output[i];
		 intptr++;
	}


	mxFree(output);

	return;
}

void Usage(void)
{
	mexErrMsgTxt("wrong usage\n");
}


