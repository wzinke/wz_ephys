#include<math.h> 
void MARgetaic(double tildA[],int nchns,int Order,int sampleno,
		     double aic[])
{
  int i,j,k;
  double det;
  for(i=0;i<Order;i++){
    j=i*(i+1)*nchns*nchns/2;
    det=1.;
    for(k=0;k<nchns;k++)
      det *= *(tildA+j+k*(nchns+1));
	/* aic[i]=-2.*((double)sampleno)*log(det)+2.*nchns*nchns*i;*/
	      aic[i]=-1.*log(det)+2.*nchns*nchns*i/((double)sampleno);
  }
  return;
}



