#include "EEGmat.h"
#include "EEGdef.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

/* 
   void MARgetpfpb(double *x[],int m,int n[],int k,int order,
		     double *Pf, double *Pb, double *Pfb)

   *x[i] is i-th channel data,
   m is number of channel.
   n[j] is the number of data in j-th segment, k is the number of 
   segments. More precisely, x[i][0:n[0]-1] input data of i channel 
   for first trial, x[i][n[0]:n[1]-1] input the second trial, and 
   x[i][n[1]:n[2]-1] input the third, and so on.

   The program return intermediate result of Pf[], Pb[], Pfb[] of certain 
   order which are used in MARfit();

   
*********************************/


void MARgetpfpb(double *x[],int m,int n[],int k,int order,
		     double Pf[],double Pb[],double Pfb[],
		     double **A, double **B)
{
  double *tempmat1;/*temp matrix*/
  double *Efn,*Ebn,*tempvec1,*tempx; /* vector */
  int i1,i2,i3,i4,i5;
  double i_n;
  tempmat1=init_mat(m);
  Efn=init_vec(m);
  Ebn=init_vec(m);
  tempvec1=init_vec(m);
  tempx=init_vec(m);

  i4=0;
  zeromat(Pf,m);
  zeromat(Pb,m);
  zeromat(Pfb,m);  
  if(order==0){
    for(i1=0;i1<k;i1++){
      for(i2=1;i2<n[i1];i2++){
	for(i3=0;i3<m;i3++)
	  tempx[i3]=x[i3][i2+i4];
	vecmvech(tempx,tempx,tempmat1,m);
	arrayp(tempmat1,Pf,Pf,m*m);
	for(i3=0;i3<m;i3++)
	  tempx[i3]=x[i3][i2+i4-1];
	vecmvech(tempx,tempx,tempmat1,m);
	arrayp(tempmat1,Pb,Pb,m*m);
      }
      i4 += n[i1];
    }
  }
  else{
    for(i1=0;i1<k;i1++){
      for(i2=order;i2<n[i1];i2++){
	zerovec(Efn,m);
	zerovec(Ebn,m);
	for(i3=0;i3<order;i3++){
	  for(i5=0;i5<m;i5++)
	    tempx[i5]=x[i5][i2+i4-i3];
	  matmvec(A[i3],tempx,tempvec1,m);
	  arrayp(tempvec1,Efn,Efn,m);
	  for(i5=0;i5<m;i5++)
	    tempx[i5]=x[i5][i2+i4-i3-1];
	  matmvec(B[order-i3-1],tempx,tempvec1,m);
	  arrayp(tempvec1,Ebn,Ebn,m);
	}
	vecmvech(Efn,Efn,tempmat1,m);
	arrayp(tempmat1,Pf,Pf,m*m);
	vecmvech(Ebn,Ebn,tempmat1,m);
	arrayp(tempmat1,Pb,Pb,m*m);
	vecmvech(Efn,Ebn,tempmat1,m);
	arrayp(tempmat1,Pfb,Pfb,m*m);
      }
      i4 += n[i1];
    }
  }

  i_n=1./((double)i4-(order+1)*k);

  matms(Pf,i_n,Pf,m);  /* normalize */
  matms(Pfb,i_n,Pfb,m);  /* normalize */
  matms(Pb,i_n,Pb,m);  /* normalize */
  free(tempmat1);
  free(Efn);
  free(Ebn);
  free(tempvec1);
  free(tempx);
  return;
  
}


