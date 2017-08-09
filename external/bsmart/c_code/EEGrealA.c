/* This is a simple transformation of tildA to realA and Ve. 
void EEGrealA(double *tildA, double **realA, double *Ve,int m, int order)

         Input: double *tildA, output of MARfit
	        int m,         channel number
                int order,     order of model you chosen
         Output:double **realA,realA matrix for furthur using in 
                               other subroutins
                double *Ve     theoretical risidual spectral matrix.

**************/

#include "EEGmat.h"
#include "EEGdef.h"

/* The interface pipe the tildA to realA, m is chanel number */ 
void EEGrealA(double *tildA, double **realA, double *Ve,int m, int order)
{
  double *A1[MAXORDER];
  int i,j;
  j=order*(order+1)*m*m/2;
  for(i=0;i<=order;i++){
    A1[i]=tildA+j;
    j += m*m;
  }
  invmat(A1[0],realA[0],m);
  matmmatH(realA[0],Ve,m);
  for(i=1;i<=order;i++)
    matmmat(realA[0],A1[i],realA[i],m);
  for(i=0;i<m*m;i++)realA[0][i]=0.;
  for(i=0;i<m;i++)realA[0][i*(m+1)]=1.;
  return;
}
