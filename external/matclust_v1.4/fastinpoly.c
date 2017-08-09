

#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>



void Usage(void) ;
double IsInPolygon(double xpoint,double ypoint,double *xvert,double *yvert,int nxvert);

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
	int               i;
	double             xmax, xmin, ymax, ymin;
        double            *xpoints;
	double            *ypoints;
	double            *xvert;
	double            *yvert;
	double			  *output;
	int				  nxpoints;
	int				  nypoints;
	int				  nxvert;
	int				  nyvert;


    /* Check numbers of arguments */
    if (!(nrhs == 4) || !(nlhs == 1)){
        Usage();
    }
   





	xpoints = mxGetPr(prhs[0]);
	nxpoints = mxGetM(prhs[0]) * mxGetN(prhs[0]);

	ypoints = mxGetPr(prhs[1]);
	nypoints = mxGetM(prhs[1]) * mxGetN(prhs[1]);

	xvert = mxGetPr(prhs[2]);
	nxvert = mxGetM(prhs[2]) * mxGetN(prhs[2]);

	yvert = mxGetPr(prhs[3]);
	nyvert = mxGetM(prhs[3]) * mxGetN(prhs[3]);

	output = (double *) mxCalloc(nxpoints, sizeof(double));

        xmax = 0;
        ymax = 0;
        xmin = 0;
        ymin = 0;
        
        for (i = 0; i < nxvert; i++) {
		
                if (xvert[i] > xmax) {
                   xmax = xvert[i];
                }
                if (yvert[i] > ymax) {
                   ymax = yvert[i];
                }
                if (xvert[i] < xmin) {
                   xmin = xvert[i];
                }
                if (yvert[i] < ymin) {
                   ymin = yvert[i];
                }                
               
         }
        
        
	for (i = 0; i < nxpoints; i++) {
		if ((xmax >= xpoints[i] >= xmin) && (ymax >= ypoints[i] >= ymin)) {
                   output[i] = IsInPolygon(xpoints[i],ypoints[i],xvert,yvert,nxvert);
                }
                else {
                   output[i] = 0;
                }
        }

	plhs[0] = mxCreateDoubleMatrix(nxpoints, 1,  mxREAL);

	dblptr = mxGetPr(plhs[0]);
	/* copy A.thetahat into the output array */
	for (i = 0; i < nxpoints; i++) {
		 *dblptr = output[i];
		 dblptr++;
	}


	mxFree(output);

	return;
}

void Usage(void)
{
	mexErrMsgTxt("wrong usage usage\n");
}


double IsInPolygon(double wx,double wy,double *xvert,double *yvert,int ncoords)
{
int	i;
int	pcross;
double	FY,bx;
    
    if(xvert == NULL) return(0);
    /*
    ** look for odd number of intersections with poly segments
    */
    pcross = 0;
    /*
    ** extend a horizontal line along the positive x axis and
    ** look at intersections
    */
    for(i=0;i<ncoords-1;i++){
	/*
	** only examine segments whose endpoint y coords
	** bracket those of the test coord
	*/
	if((wy > yvert[i] && wy <= yvert[i+1]) ||
	(wy < yvert[i] && wy >= yvert[i+1])){
	    /*
	    ** count those which are on the positive x side 
	    ** by computing the intercept.
	    ** find the x value of the line between the two fcoords
	    ** which is at wy
	    */
/*
	    (yvert[i] - wy)/(xvert[i] - bx) = 
		(yvert[i+1] - wy)/(xvert[i+1] - bx);

	    (f2x - bx)/(f1x - bx) = (f2y - wy)/(f1y - wy);
	    (f2x - bx) = (f1x - bx)*(f2y - wy)/(f1y - wy);
	    bx((f2y-wy)/(f1y-wy) -1)= f1x*(f2y-wy)/(f1y-wy) - f2x;
	    bx = (f1x*(f2y-wy)/(f1y-wy)-f2x)/((f2y-wy)/(f1y-wy)-1)

	    FY = (f2y-wy)/(f1y-wy);
	    bx = (f1x*FY - f2x)/(FY-1);
*/
	    FY = (yvert[i+1]-wy)/(yvert[i]-wy);
	    if(FY == 1){
		pcross++;
	    } else {
		bx = (xvert[i]*FY - xvert[i+1])/(FY-1);

		if(bx >= wx){
		    pcross++;
		}
	    }
	}
    }
    if(i == ncoords-1){
	/*
	** compute the final point which closes the polygon
	*/
	if((wy > yvert[i] && wy <= yvert[0]) ||
	(wy < yvert[i] && wy >= yvert[0])){
	    FY = (yvert[0]-wy)/(yvert[i]-wy);
	    if(FY == 1){
		pcross++;
	    } else {
		bx = (xvert[i]*FY - xvert[0])/(FY-1);
		if(bx > wx){
		    pcross++;
		}
	    }
	}
    }
    /*
    ** now look for an odd number of crossings
    */
    if(pcross > 0 && (pcross%2 == 1)){
	return(1);
    } else {
	return(0);
    }
}
