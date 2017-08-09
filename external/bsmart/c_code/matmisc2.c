#include <math.h>
#include <stdlib.h>
#define TINY 1.0e-20;
double *init_vec();
double *init_mat();
#include<stdio.h>
void EEGerror(char *s)
{
  fprintf(stderr,"%s\n",s);
  exit(1);
}

void ludcmp(double *a, int n, int *indx, double *d)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;

	vv=init_vec(n);
	*d=1.0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i*n+j])) > big) big=temp;
		if (big == 0.0) EEGerror("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[i*n+j];
			for (k=0;k<i;k++) sum -= a[i*n+k]*a[k*n+j];
			a[i*n+j]=sum;
		}
		big=0.0;
		for (i=j;i<n;i++) {
			sum=a[i*n+j];
			for (k=0;k<j;k++)
				sum -= a[i*n+k]*a[k*n+j];
			a[i*n+j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<n;k++) {
				dum=a[imax*n+k];
				a[imax*n+k]=a[j*n+k];
				a[j*n+k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j*n+j] == 0.0) a[j*n+j]=TINY;
		if (j != n-1) {
			dum=1.0/(a[j*n+j]);
			for (i=j+1;i<n;i++) a[i*n+j] *= dum;
		}
	}
	free(vv);
}
#undef TINY

void lubksb(double *a, int n, int *indx, double b[])
{
	int i,ii=-1,ip,j;
	double sum;

	for (i=0;i<n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii!=-1)
			for (j=ii;j<=i-1;j++) sum -= a[i*n+j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[i*n+j]*b[j];
		b[i]=sum/a[i*n+i];
	}
}

/* inverse general matrix */
void invmat(double *A, double *B, int m)
{
  double d,*col,*temp;
  int i,j,*indx;
  col=init_vec(m);
  temp=init_mat(m);
  if((indx=malloc(m*sizeof(int)))==NULL)
    EEGerror("invmat---memory allocation problem\n");
  arraycopy(A,temp,m*m);
  ludcmp(temp,m,indx,&d);
  for(j=0;j<m;j++){
    for(i=0;i<m;i++)col[i]=0.0;
    col[j]=1.0;
    lubksb(temp,m,indx,col);
    for(i=0;i<m;i++)B[i*m+j]=col[i];
  }
  free(col);
  free(temp);
  free(indx);
  return;
}

/* determinant of matrix */
void detmat(double *A, double *x, int m)
{
  double *temp,d;
  int i,*indx;
  temp=init_mat(m);
  if((indx=malloc(m*sizeof(int)))==NULL)
    EEGerror("detmat--memory allocation problem\n");
  arraycopy(A,temp,m*m);
  ludcmp(temp,m,indx,&d);
  for(i=0;i<m;i++) d *= temp[i*m+i];
  *x=d;
  free(temp);
  free(indx);
  return;
}

/* This subroutine calculate the square root matrix B of matrix A
the size of matrix is n*n. the matrix A is symmetric, the matrix B is 
lower triangle, and A=B*B^H. For avoiding error on A[][], I use one 
dimensional array to represent matrix, the order A[i][j]=A[m*i+j]*/

void sqrmat(double *A, double *B, int m)
{
  int i,j,k;
  double x;
  for(i=0;i<m-1;i++)
    for(j=i+1;j<m;j++)
      B[m*i+j]=0.;
  for(i=0;i<m;i++){
    x=0.;                /* B[i][i]=sqrt(A[i][i]-/sum_{j=0}^{i-1}B[i][j]) */
    for(j=0;j<i;j++)
      x += B[i*m+j]*B[i*m+j];
    x=A[m*i+i]-x;
    if(x<0.)EEGerror("sqrmat--matrix A is not positive definite");
    B[m*i+i]=sqrt(x);
    for(j=i+1;j<m;j++){ /* B[j][i]=(A[j][i]-\sum_{k=0}^{i-1}B[i][k]*B[j][k])
			           /B[i][i]  */
      x=0.;
      for(k=0;k<i;k++)
	x += B[i*m+k]*B[j*m+k];
      x=A[j*m+i]-x;
      B[j*m+i]=x/B[m*i+i];
    }
  }
  return;
}

/* The subroutine reverse a lower triangle matrix A to matrix B
          B=A^{-1}, matrix m*m, one dimensional array represent matrix*/

void revltmat(double *A, double *B, int m)
{
  int i,j,k;
  double x;
  for(i=0;i<m-1;i++)
    for(j=i+1;j<m;j++)
      B[m*i+j]=0.;
  for(i=0;i<m;i++)
    B[m*i+i]=1./A[m*i+i];    /* B[i][i]=1./A[i][i] */
  for(i=1;i<m;i++){/*B[j][j-i]=-(\sum_{k=j-i}^{j-1}A[j][k]B[k][j-i])/A[j][j]*/
    for(j=i;j<m;j++){
      x=0.;
      for(k=j-i;k<j;k++)
	x -= A[j*m+k]*B[k*m+j-i];
      B[j*m+j-i]= x/A[j*m+j];
    }
  }
  return;
}

void zerovec(double *x, int n)
{
  int i;
  for(i=0;i<n;i++)
    x[i]=0.0;
}
