#include <math.h>
#define MAX(a,b) ((a)>(b)? (a) : (b))
#define MIN(a,b) ((a)<(b)? (a) : (b))


void minmax(long nmax, int *x, int *xmin, int *xmax)
{
   long i;

   for(*xmin=x[0], *xmax=x[0], i=0; i<nmax; i++){
       *xmax=MAX(*xmax,x[i]);
       *xmin=MIN(*xmin,x[i]);
   }
}
/****
void minmaxdouble(long nmax, double *x, double *xmin, double *xmax)
{
   long i;

   for(*xmin=x[0], *xmax=x[0], i=0; i<nmax; i++){
       *xmax=MAX(*xmax,x[i]);
       *xmin=MIN(*xmin,x[i]);
   }
}

*****/



double cor(double *sigarray1, double mean1, double ss1, double *sigarray2,
double mean2, double ss2, int npoints, int nlags, double *corfun)
{
  int lag;
  double cov, getcov(), sqrt();

  for(lag = 1; lag <= nlags; lag++)
    {
      cov = getcov(sigarray1,sigarray2,mean1,mean2,npoints,lag);
      corfun[lag-1] = cov / sqrt(ss1 * ss2);
    }
  
  for(lag = 1; lag <= nlags; lag++)
    {
      cov = getcov(sigarray2,sigarray1,mean2,mean1,npoints,lag);
      corfun[lag+nlags-1] = cov / sqrt(ss1 * ss2);
    }
}

double getmean(double *array, int npoints)
{
  int i;
  double mean = 0.0;

  for (i = 0; i < npoints; i++)
    {
      mean += array[i]; 
    }

  mean /= (double)npoints;

  return(mean);
}

double getss(double *array, double mean, int npoints)
{
  int i;
  double ss = 0.0;

  for (i = 0; i < npoints; i++)
    {
      ss += (array[i] - mean) * (array[i] - mean);
    }

  return(ss);
}

double getcov(double *array1, double *array2, double mean1, double mean2, int npoints, int lag)
{
  int i;
  double covar = 0.0;

  for(i=0; i < npoints-lag; i++)
    {
      covar += (array1[i] - mean1) * (array2[i+lag] - mean2);
    }

  return(covar);

}




