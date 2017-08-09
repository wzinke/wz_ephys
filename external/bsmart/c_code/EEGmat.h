#ifndef EEGMAT_H
#define EEGMAT_H

/* matmisc.c */
double* init_vec(int);
double* init_mat(int);
void vecmvech(double *, double *, double *, int );
void matmvec(double *, double *, double *, int );
void arraym(double *, double *, double *, int );
void matH(double *, double *, int );
void arraycopy(double *, double *, int );
void matmmat(double *, double *, double *, int );
void ltmatmmat(double *, double *, double *, int );
void matmutmat(double *, double *, double *, int );
void matmmatH(double *, double *, int );
void arrayp(double *, double *, double *, int );
void zeromat(double *,int );
void Immat(double *,double *, int);
void matms(double* ,double, double*, int);
/* matmisc2.c */
void EEGerror(char *);
void invmat(double *, double *, int );
void detmat(double *, double *, int );
void sqrmat(double *, double *, int );
void revltmat(double *, double *, int );
void zerovec(double *,int );

/*matmisc5.c*/
void minmax(long , int *, int *, int *);
double cor(double *, double, double, double *,
		   double, double, int, int, double *);
double getcov(double *, double *, double, double, int, int);
double getmean(double *, int);
double getss(double *, double, int);

/*matmisc3.c*/
typedef struct {
  double x;  /* Re(C) */
  double y;  /* Im(C) */
} Complex;

double sCabs(Complex);
double Cabs(Complex);
Complex Cadd(Complex , Complex );
Complex complex(double ,double );
void invCmat(Complex *A, Complex *B, int m);
Complex Cmul(Complex, Complex);
void CmatmCmatH(Complex *A,Complex *B,int m);
void CmatHmCmat(Complex *A,Complex *B,int m);
void detCmat(Complex *A,Complex *x,int m);
Complex Conj(Complex);
void CmatmCmat(Complex *A,Complex *B,Complex *C,int m);
void ConjCmat(Complex *A,Complex *B,int m);
void Carrayp(Complex *A, Complex *B, Complex *C, int m);
#endif
