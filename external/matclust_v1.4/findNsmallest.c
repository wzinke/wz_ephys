

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

      double 			  *dblptr;
      double      *values;
    double         currentLargest;
    long int        currentLargestInd;
      
      double      *numValuesToGetPr;
      long int      numValuesToGet;
      
      double			  *output1;
    double              *output2;
      
    long int				  nvalues,i, j, currentCheckVal;
    bool                 filled = false;
    
    /* Check numbers of arguments */
    if (!((nrhs == 2)) || !((nlhs == 1)||(nlhs == 0)||(nlhs == 2)) ) {
         Usage();
    }
   

    values = mxGetData(prhs[0]);
    nvalues = mxGetM(prhs[0]) * mxGetN(prhs[0]);
    
    if ((mxGetM(prhs[1]) * mxGetN(prhs[1])) != 1) {
        mexErrMsgTxt("numValues must be a scalar \n");
    }

    numValuesToGetPr = mxGetData(prhs[1]);
    numValuesToGet = numValuesToGetPr[0];
    if (numValuesToGet < 1) {
        mexErrMsgTxt("numValues must be a positive integer \n");

    }
    if (numValuesToGet > nvalues) {
        mexErrMsgTxt("numValues must be less than the length of values \n");
        
    }

      
                     
    output1 = (double *) mxCalloc(numValuesToGet, sizeof(double));
    output2 = (double *) mxCalloc(numValuesToGet, sizeof(double));
    
    currentCheckVal = 0;
    currentLargest = values[0]-1;

    for (i = 0; i < nvalues; i++) {
        //mexPrintf(i); 
        //mexEvalString("here"); 
        if (filled) {
            
            if (values[i] < currentLargest) {
                output1[currentLargestInd] = values[i];
                output2[currentLargestInd] = i+1;
            
            
                //currentLargest = values[i];
                currentLargest = output1[0];
                currentLargestInd = 0;
                
                for (j=0; j < numValuesToGet; j++) {
                    if (output1[j] > currentLargest)  {
                        currentLargest = output1[j];
                        currentLargestInd = j;
                    }
                }
            }
    
            
        } else {
            output1[currentCheckVal] = values[i];
            output2[currentCheckVal] = i+1;
            if (values[i] > currentLargest) {
                currentLargest = values[i];
                currentLargestInd = currentCheckVal;
            }
            currentCheckVal++;
            if (currentCheckVal >= numValuesToGet) {
                filled = true;
            }
            
        }
    }
             
             
                   
                                                      
    plhs[0] = mxCreateDoubleMatrix(numValuesToGet,1, mxREAL);
    dblptr = mxGetPr(plhs[0]);
    for (i = 0; i < numValuesToGet; i++) {
            *dblptr = output1[i];
            dblptr++;
    }
    
    if (nlhs > 1) {
        plhs[1] = mxCreateDoubleMatrix(numValuesToGet,1, mxREAL);
        dblptr = mxGetPr(plhs[1]);
        for (i = 0; i < numValuesToGet; i++) {
            *dblptr = output2[i];
            dblptr++;
        }
    }
  
    
        
    return;
}

void Usage(void)
{
   mexErrMsgTxt("Usage: 2 inputs and 1 or 2 outputs\n");
}


