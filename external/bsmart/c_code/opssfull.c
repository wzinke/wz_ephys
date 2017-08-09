/*   Usage: opss datafile ARcoeff Noisefile   
   for one given widow, parameters need modify, see
   EEGdef.h,             */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "EEGdef.h"
#include "EEGmat.h"


main(int argc, char**argv)
{ 
  int NPTS =0;
  int NCHN =0;
  int NTRLS =0;
  int WIN=0;
  int MODORDER=0;
if (argc==9)
   {
     NCHN=atoi(argv[5]);
     NTRLS=atoi(argv[6]);
     NPTS=atoi(argv[7]);
     MODORDER=atoi(argv[8]);
     WIN=NPTS*NTRLS;
   }
else if ( argc!=5 && argc!=9)
	{
	  fprintf(stderr,"opssfull datfile A Ve AIC Nchannels Ntrails Npoints MODORDER\n");
	  exit(1);
	}
else { 
  double chan[1],trai[1],poin[1],order[1];
  FILE *chanfp,*traifp,*poinfp,*windfp,*orderfp;
  if((chanfp=fopen("channel","r"))==NULL)   
              printf("The file'channel' was not opened\n");   
  else   
    fscanf(chanfp,"%lf",&chan[0]);  
    NCHN=(int)chan[0];
    fclose(chanfp);
  if((traifp=fopen("trail","r"))==NULL)   
              printf("The file'trail' was not opened\n");   
  else   
    fscanf(traifp,"%lf",&trai[0]);  
    NTRLS=(int)trai[0];
    fclose(traifp);
  if((poinfp=fopen("points","r"))==NULL)   
              printf("The file'points' was not opened\n");   
  else   
    fscanf(poinfp,"%lf",&poin[0]);  
    NPTS=(int)poin[0];
    fclose(poinfp);
  if((orderfp=fopen("order","r"))==NULL)   
              printf("The file'order' was not opened\n");   
  else   
    fscanf(orderfp,"%lf",&order[0]);  
    MODORDER=(int)order[0];
    fclose(orderfp);
  WIN=NPTS*NTRLS;
}
  FILE *inpt, *cpt, *ppt, *fp;

  char outfil[100], cohfil[100], phsfil[100];
  char chin1[3], chin2[3];

  int ch1, ch2, index, getindex();

  double *A[MAXORDER],*Ve,*tildA;
  double *comat;
  double freq, phase;
  double **x;

  Complex *H,*S,*P;

  float **dat; 
  int *n;  /* n[j] is the number of data in j-th segment, 120 */  
  int i, j, trial, chn;
  int idx=0;
 
  double aic[20];



  for(i=0;i<MAXORDER;i++){
    if((A[i]=malloc(NCHN*NCHN*sizeof(double)))==NULL)
      EEGerror("main---memory allocation error\n");
  }

  if((Ve=malloc(NCHN*NCHN*sizeof(double)))==NULL)
    EEGerror("main---memory allocation error\n");
  if((comat=malloc(NCHN*NCHN*sizeof(double)))==NULL)
    EEGerror("main---memory allocation error\n");
  if((H=malloc(NCHN*NCHN*sizeof(Complex)))==NULL)
    EEGerror("main---memory allocation error\n");
  if((S=malloc(NCHN*NCHN*sizeof(Complex)))==NULL)
    EEGerror("main---memory allocation error\n");
  if((P=malloc(NCHN*NCHN*sizeof(Complex)))==NULL)
    EEGerror("main---memory allocation error\n");
  if((tildA=malloc(MAXORDER*(MAXORDER+1)*NCHN*NCHN*sizeof(double)/2))==NULL)
    EEGerror("main---memory allocation error\n");



  /* allocation of memory for dat */
  dat = malloc(NCHN*sizeof(float*));
  for( i = 0; i < NCHN; i++)
    dat[i] = malloc(NPTS*NTRLS*sizeof(float));

  x = malloc(NCHN*sizeof(double*));
  for( i = 0; i < NCHN; i++)
    x[i] = malloc(NPTS*NTRLS*sizeof(double));

  /* required by MARfit */ 	
   n=malloc(NTRLS*sizeof(int));
  for( i = 0; i < NTRLS; i++) n[i] = NPTS;
  
   inpt = fopen(argv[1],"rb"); 
/*  inpt = fopen("gecesar4abcde.dat","rb"); */

  for( trial = 0; trial < NTRLS; trial++){
    for( chn = 0; chn < NCHN; chn++){

	  if( fread(&dat[chn][idx],sizeof(float),NPTS,inpt) !=NPTS) {
		if(feof(inpt)) {printf("premature end of file\n");exit(-1);}
		else {printf("file read error\n");exit(-1);}
		
	  }
   
	}
   idx+=NPTS;
    /* printf("ok = %d\n", idx); */
  }  /* end of trial */	


/* convert data format from float to double  */
  for (i=0; i < NCHN; i++)
	for(j=0; j < NPTS*NTRLS; j++)
	  x[i][j] = dat[i][j];
 /*printf("y2k\n");*/

    MARfit(x,NCHN,n,NTRLS,MODORDER,tildA); 
  /*  MARfit(x,NCHN,n,NTRLS,6,tildA);  */

/* AIC to determine the model order */
  j=0;
 for(i=0;i<NTRLS;i++) j += n[i];
 MARgetaic(tildA,NCHN,20,j,aic);
   fp = fopen(argv[4],"w"); 
   for(i=0;i<20;i++)
	 /* fprintf(stderr,"%d  %g\n",i,aic[i]);  */
	 fprintf(fp,"%d  %g\n",i,aic[i]);
   fclose(fp);
//printf("ok\n");

  EEGrealA(tildA,A,Ve,NCHN,MODORDER);
//printf("oook\n");

   fp = fopen(argv[2],"w"); 
   for ( i=0; i < MODORDER+1; i++)   /* org=6,if 7, then new =0 */
	 for ( j=0; j < NCHN*NCHN; j++)
	   fprintf(fp,"%.3g ",A[i][j]);
   fclose(fp);
   
   fp = fopen(argv[3],"w"); 
   for ( i=0; i < NCHN*NCHN; i++)
	 fprintf(fp,"%.3g ",Ve[i]);
   fclose(fp);
   


 /* for(freq=0.;freq<1.;freq += 0.01){
	MAR_h(A,H,NCHN,5,freq);

	  save transfer function H 
   fp = fopen(argv[4],"a+"); 
   for ( i=0; i < NCHN*NCHN; i++)
         fprintf(fp,"%.3g %.3g\n",H[i].x, H[i].y);
   fclose(fp);
   

	MAR_s(H,Ve,S,NCHN);
	MAR_comat(S,comat,NCHN);
	MAR_sphase(S,P,NCHN);
	
	for (ch1 = 0; ch1 < NCHN; ch1++)
	  {
	    for (ch2 = 0; ch2 < NCHN; ch2++)
	      {
		if(ch2 > ch1)
		  {
		    index = getindex(ch1,ch2,NCHN);

		    itoa(chin1,ch1+1);
		    itoa(chin2,ch2+1);
		    strcpy(outfil,chin1);
		    strcat(outfil,"_");
		    strcat(outfil,chin2);

		    strcpy(cohfil,outfil);
		    strcpy(phsfil,outfil);

		    strcat(cohfil,".coh");
		    strcat(phsfil,".phs");

		    cpt=fopen(cohfil,"a+");
		    fprintf(cpt,"%.3g %.3g\n",freq,comat[index]);
		    fclose(cpt);

		    phase = atan2(P[index].y,P[index].x);
		    ppt=fopen(phsfil,"a+");
		    fprintf(ppt,"%.3g %.3g\n",freq,phase);
		    fclose(ppt);
		  }
	      }
	  }
      }  end of freq loop */


  free(tildA);
  free(Ve);
  free(comat);
  free(H);
  free(S);
  free(P);
  for(i=0;i<MAXORDER;i++)free(A[i]);
  free(n);

  for (i = 0; i < NCHN; i++)
    free(dat[i]);
  free(dat);

  for (i = 0; i < NCHN; i++)
    free(x[i]);
  free(x);

  exit(0);


   
}


/*int getindex(int ch1,int ch2)
{
  int index;

  index = (ch1 * NCHN) + ch2;
  return(index);
}*/





