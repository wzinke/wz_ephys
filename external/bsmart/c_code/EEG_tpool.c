/* #include<math.h>  */
#include<stdlib.h>
void EEG_tpool(double **xi, int *ni, int ki, int nchns,
	       double **xo, int *no, int ko)
{
  int i,j,k,pos,a_pos_in,a_pos_out;
  a_pos_out=0;
  /* printf("in and out number = %d  %d\n", ki, ko); */
  for(i=0;i<ko;i++){
    pos = (int) ( ((double)ki)*drand48() );
	/* printf("dddd = %d\n", pos); */
    while(pos>=ki || pos<0)         /* make sure 0 <= pos < ki */    
      pos = (int) ( ((double)ki)*drand48() );
    no[i]=ni[pos];
    a_pos_in=0;
    for(j=0;j<pos;j++)
      a_pos_in += ni[j];
    for(j=0;j<nchns;j++){
      for(k=0;k<no[i];k++){
	xo[j][k+a_pos_out] = xi[j][k+a_pos_in];
      }
    }
    a_pos_out += no[i];
  }
  return;
}
