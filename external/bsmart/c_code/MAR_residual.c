#include "EEGmat.h"
void MAR_residual(double **xin, double **rout, double **A,
	     int nchns, int order, int xinlength)
{
  int i,j,k;
  double *temp,*y1,*y2;
  temp=init_vec(nchns);
  y1=init_vec(nchns);
  y2=init_vec(nchns);
  for(i=0;i<xinlength-order;i++){
    for(j=0;j<nchns;j++)y1[j]=0.;
    for(j=0;j<=order;j++){
      for(k=0;k<nchns;k++)
	y2[k]=xin[k][i+j];
      matmvec(A[j],y2,temp,nchns);
      arrayp(temp,y1,y1,nchns);
    }
    for(j=0;j<nchns;j++)
      rout[j][i]=y1[j];
  }
  free(temp);
  free(y1);
  free(y2);
  return;
}
      
      
