/* matrix calculation subroutine */

#include<stdlib.h>
void EEGerror();

/* 1. initial vector */

double* init_vec(int n)
{
  double *x;
  int i;
  
  if((x=malloc(n*sizeof(double)))==NULL)
    EEGerror("init_vec--memory allocation problem\n");
  for(i=0;i<n;i++)
    x[i]=0.;
  return(x);
}

/* 2. initial matrix */

double* init_mat(int m)
{
  double *x;
  int i;
  
  if((x=malloc(m*m*sizeof(double)))==NULL)
    EEGerror("init_mat--memory allocation problem\n");
  for(i=0;i<m*m;i++)
    x[i]=0.;
  return(x);
}

/* 3. vec * vec^H = mat */

void vecmvech(double *a, double *b, double *C, int m)
{
  int i,j;
  for(i=0;i<m;i++)
    for(j=0;j<m;j++)
      C[i*m+j]=a[i]*b[j];
  return;
}

/* 4. mat * vec = vec */

void matmvec(double *A, double *b, double *c, int m)
{
  int i,j;
  for(i=0;i<m;i++){
    c[i]=A[i*m]*b[0];
    for(j=1;j<m;j++)
      c[i] += A[i*m+j]*b[j];
  }
  return;
}

/* 5. array - array = array */

void arraym(double *a, double *b, double *c, int n)
{
  int i;
  for(i=0;i<n;i++)
    c[i]=a[i]-b[i];
  return;
}

/* 6. Transpose of mat */

void matH(double *A, double *B, int m)
{
  int i,j;
  for(i=0;i<m;i++)
    for(j=0;j<m;j++)
      B[i*m+j]=A[j*m+i];
  return;
}

/* 7. array = array */

void arraycopy(double *A, double *B, int m)
{
  int i;
  for(i=0;i<m;i++)
    B[i]=A[i];
  return;
}

/*  8. mat * mat = mat */

void matmmat(double *A, double *B, double *C, int m)
{
  int i,j,k,n;
  for(i=0;i<m;i++)
    for(j=0;j<m;j++){
      n=i*m+j;
      C[n]=0.;
      for(k=0;k<m;k++)
	C[n] += A[i*m+k]*B[k*m+j];
    }
  return;
}

/* 9. ltmat * mat = mat */

void ltmatmmat(double *A, double *B, double *C, int m)
{
  int i,j,k,n;
  for(i=0;i<m;i++)
    for(j=0;j<m;j++){
      n=i*m+j;
      C[n]=0.;
      for(k=0;k<=i;k++)
	C[n] += A[i*m+k]*B[k*m+j];
    }
  return;
}

/* 10. mat * utmat = mat */

void matmutmat(double *A, double *B, double *C, int m)
{
  int i,j,k,n;
  for(i=0;i<m;i++)
    for(j=0;j<m;j++){
      n=i*m+j;
      C[n]=0.;
      for(k=0;k<=j;k++)
	C[n] += A[i*m+k]*B[k*m+j];
    }
  return;
}

/* 11 mat * matH = mat */

void matmmatH(double *A, double *B, int m)
{
  int i,j,k,n;
  for(i=0;i<m;i++)
    for(j=i;j<m;j++){
      n=i*m+j;
      B[n]=0.;
      for(k=0;k<m;k++)
	B[n] += A[i*m+k]*A[j*m+k];
    }
  for(i=1;i<m;i++)
    for(j=0;j<i;j++)
      B[i*m+j]=B[j*m+i];
  return;
}

/*  12. array + array = array */

void arrayp(double *a, double *b, double *c, int n)
{
  int i;
  for(i=0;i<n;i++)
    c[i]=a[i]+b[i];
  return;
}

/* 13. matrix=0 */

void zeromat(double *A, int m)
{
  int i,j;
  j=m*m;
  for(i=0;i<j;i++)
    A[i]=0.0;
}

/* 14. B = I - A  */

void Immat(double *A, double *B, int m)
{
  int i;
  for(i=0;i<m*m;i++)
    B[i]=-A[i];
  for(i=0;i<m;i++)
    B[i*(m+1)] += 1.;
}

/* 15. matrix A * scalar x = matrix B*/

void matms(double *A, double x,double *B, int m)
{
  int i,j;
  j=m*m;
  for(i=0;i<j;i++)
    B[i] = A[i]*x;
  return;
}
