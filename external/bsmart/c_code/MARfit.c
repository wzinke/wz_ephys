/*** Content *******

This program is written for realizing the normalized LWR algorithm for 
Multivariate Autoregressive model (reference: S. Haykin and S.Kesler: 
"Prediction-error filtering and Maximum-Entropy Spetral Estimation" 
in "Nonlinear Methods of Spectral Analysis" ed. S. Haykin, pp 9-72(at pp58),
Springer-Verlag, New York 1983, and Martin Morf, Augusto Vieira, 
Daniel T. L. Lee, Thomas Kailath, "Recursive Multichannel Maximum Entropy 
Spectral Estimation", IEEE Trans. GE-16, No. 2, 85-94(1978). ) 

I can't realize FPE criteria, so I only use AIC criteria, which is 
determinant of residual noise. 

*/

#include "EEGmat.h"
#include "EEGdef.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

/* 
   void MARfit(double *x[],int m,int n[],int k,int Order,
		     double tildeA[])

   Solving MAR model  for p=0 to Order(up to MAXORDER which is 50 that defined
            in EEGdef.h), 
   
   *x[i] is i-th channel data,
   m is number of channel.
   n[j] is the number of data in j-th segment, k is the number of 
   segments. More precisely, x[i][0:n[0]-1] input data of i channel 
   for first trial, x[i][n[0]:n[1]-1] input the second trial, and 
   x[i][n[1]:n[2]-1] input the third, and so on.

   The program return tildeA[Order*(Order+1)*m*m/2]
   In form of tildA[ord*(ord+1)*m*m/2+t*m*m+i*m+j] it outputs {ij} element 
   of tildA_t in "ord" order model.

   Dependence:
        MARgetpfpb();
   
*********************************/


void MARfit(double *x[],int m,int n[],int k,int Order,
		     double tildeA[])
{
  double *A[MAXORDER],*B[MAXORDER],*Pf,*Pb,*Pfb,*rho,*rhoH,*P,*Q; /* matrix */
  double *tempmat1,*tempmat2,*tempmat3,*tempmat4,*tempmat5[MAXORDER];/*temp matrix*/
  double *tempA[50];  /* temp pointer */
  int i1,i2,order;
  int tAcount;
  Order = Order + 1;
  for(i1=0;i1<=Order;i1++){
    A[i1]=init_mat(m);
    B[i1]=init_mat(m);
    tempmat5[i1]=init_mat(m);
  }
  Pf=init_mat(m);
  Pb=init_mat(m);
  Pfb=init_mat(m);
  P=init_mat(m);
  Q=init_mat(m);
  rho=init_mat(m);
  rhoH=init_mat(m);
  tempmat1=init_mat(m);
  tempmat2=init_mat(m);
  tempmat3=init_mat(m);
  tempmat4=init_mat(m);
  tAcount=0;

  MARgetpfpb(x,m,n,k,0,Pf,Pb,Pfb,A,B);

/*    A(0)   */
  sqrmat(Pf,tempmat1,m);
  revltmat(tempmat1,A[0],m);
  for(i1=0;i1<m*m;i1++){
    tildeA[tAcount]=*(A[0]+i1);
    tAcount++;
  }

/*      B(0)    */
  sqrmat(Pb,tempmat1,m);
  revltmat(tempmat1,B[0],m);

  for(order=1;order<Order;order++){

    MARgetpfpb(x,m,n,k,order,Pf,Pb,Pfb,A,B);

/* rho=Pf^(-1/2) Pfb Pb^(-H/2)  */

    sqrmat(Pf,tempmat1,m);
    revltmat(tempmat1,tempmat2,m);
    sqrmat(Pb,tempmat1,m);
    revltmat(tempmat1,tempmat3,m);
    matH(tempmat3,tempmat1,m);
    ltmatmmat(tempmat2,Pfb,tempmat3,m);
    matmutmat(tempmat3,tempmat1,rho,m);

/* P=I-rho*rhoH  Q=I-rhoH * rho    */
    matH(rho,rhoH,m);
    matmmatH(rho,tempmat1,m);
    matmmatH(rhoH,tempmat2,m);
    Immat(tempmat1,P,m);
    Immat(tempmat2,Q,m);

/*   update A[i] and B[i]  */
    for(i1=0;i1<order;i1++)
      arraycopy(A[i1],tempmat5[i1],m*m);
    sqrmat(P,tempmat3,m);
    revltmat(tempmat3,tempmat1,m);
    sqrmat(Q,tempmat3,m);
    revltmat(tempmat3,tempmat2,m);

    /* A[0]  */
    ltmatmmat(tempmat1,A[0],tempmat3,m);
    arraycopy(tempmat3,A[0],m*m);

    for(i1=1;i1<order;i1++){
      matmmat(rho,B[order-i1],tempmat3,m);
      arraym(A[i1],tempmat3,tempmat4,m*m);
      ltmatmmat(tempmat1,tempmat4,A[i1],m);
    }

    /* A[order] */
    matmmat(rho,B[0],tempmat3,m);
    matms(tempmat3,-1,tempmat3,m);
    ltmatmmat(tempmat1,tempmat3,A[order],m);

    for(i1=0;i1<=order;i1++){
      for(i2=0;i2<m*m;i2++){
        tildeA[tAcount]=*(A[i1]+i2);
        tAcount++;
      }
    }

    /* B[0]  */
    ltmatmmat(tempmat2,B[0],tempmat3,m);
    arraycopy(tempmat3,B[0],m*m);

    for(i1=1;i1<order;i1++){
      matmmat(rhoH,tempmat5[order-i1],tempmat3,m);
      arraym(B[i1],tempmat3,tempmat4,m*m);
      ltmatmmat(tempmat2,tempmat4,B[i1],m);
    }

    /* B[order] */
    matmmat(rhoH,tempmat5[0],tempmat3,m);
    matms(tempmat3,-1.,tempmat3,m);
    ltmatmmat(tempmat2,tempmat3,B[order],m);

  }
  
/* Free memory */
  for(i1=0;i1<Order;i1++){
    free(A[i1]);
    free(B[i1]);
    free(tempmat5[i1]);
  }
  free(Pf);
  free(Pb);
  free(Pfb);
  free(P);
  free(Q);
  free(rho);
  free(rhoH);
  free(tempmat1);
  free(tempmat2);
  free(tempmat3);
  free(tempmat4);
  return;
  
}


