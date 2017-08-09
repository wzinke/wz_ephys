#include <R.h>
#include <Rmath.h>
#include <math.h>

#define RECT(x) x<0.0 ? 0.0 : x 

static double cFct(double t, double alpha, double beta);
static double bFct(double t, double beta);
static double rhFct(double t, double alpha, double beta);
static double logrhFct(double t, double alpha, double beta);
static double kFct(double t, double u, double alpha, double beta);
static double logkFct(double t, double u, double alpha, double beta);

void crossTight(double *tMax, int *n, double *a, double *b, double *G) {

  int i,j;
  double h,sum,t,u;
  double k;
  
  h = *tMax / *n;
  t = h;
  u = 0.5*h;
  G[0] = RECT(rhFct(t,*a,*b)/kFct(t,u,*a,*b));
  
  for (j = 2; j <= *n; j++) {
    t = j*h;
    sum = rhFct(t,*a,*b);
    
    for (i = 1; i < j; i++) {
      u = (i-0.5)*h;
      k = kFct(t,u,*a,*b);
      sum -= k*G[i-1];
    }
    u = (j-0.5)*h;
    k = kFct(t,u,*a,*b);
    G[j-1] = RECT(sum/k);
  } 
  
  for (i = 1; i < *n; i++) {
    G[i] += G[i-1];
  }

} 

void crossTightWithB(double *tMax, int *n, double *a, double *b, double *G, double *Gu, double *Gl) {

  int i,j;
  double h,sum,t,u;
  double f,k,kj,kjp,kjm;
  double sumU,sumL;
  
  h = *tMax / *n;
  t = h;
  u = 0.5*h;
  f = rhFct(t,*a,*b);
  Gl[0] = f;
  G[0] = RECT(f/kFct(t,u,*a,*b));
  Gu[0] = f/kFct(t,0.0,*a,*b);
  
  for (j = 2; j <= *n; j++) {
    t = j*h;
    sum = rhFct(t,*a,*b);
    
    for (i = 1; i < j; i++) {
      u = (i-0.5)*h;
      k = kFct(t,u,*a,*b);
      sum -= k*G[i-1];
    }
    u = (j-0.5)*h;
    k = kFct(t,u,*a,*b);
    G[j-1] = RECT(sum/k);
  } 
  
  for (i = 1; i < *n; i++) {
    G[i] += G[i-1];
  }

  for (i = 1; i < *n; i++) {
    t = (i+1)*h;
    f = rhFct(t,*a,*b);
    sumU = f;
    sumL = f;
    
    kjm = kFct(t,0.0,*a,*b);
    kj = kFct(t,h,*a,*b);
    for (j = 0; j < i; j++) {
      if (j < i-1) {
	u = (j+1)*h;
	kjp = kFct(t,u+h,*a,*b);
      } else {
	kjp = 1.0;
      }
      sumL += (kjp-kj)*Gl[j];
      sumU += (kj-kjm)*Gu[j];
      kjm=kj;
      kj=kjp;
    }
    Gl[i] = sumL;
    Gu[i] = sumU/kjm;
  } 
  
} 


void crossTightlog(double *tMax, int *n, double *a, double *b, double *G) {

  int i,j;
  double h,sum,t,u;
  double k;
  
  h = *tMax / *n;
  t = h;
  u = 0.5*h;
  G[0] = logrhFct(t,*a,*b) - logkFct(t,u,*a,*b);;
  
  for (j = 2; j <= *n; j++) {
    t = j*h;
    sum = logrhFct(t,*a,*b);
    
    for (i = 1; i < j; i++) {
      u=(i-0.5)*h;
      k=logkFct(t,u,*a,*b);
      sum = logspace_sub(sum,k+G[i-1]);
    }
    u=(j-0.5)*h;
    k=logkFct(t,u,*a,*b);
    G[j-1] = sum-k;
  } 
  
  for (i = 1; i < *n; i++) {
    G[i] = logspace_add(G[i-1],G[i]);
  }
  
  for (i = 0; i < *n; i++) {
    G[i] = exp(G[i]);
  }

} 


void crossTightWithBlog(double *tMax, int *n, double *a, double *b, double *G, double *Gu, double *Gl) {

  int i,j;
  double h,sum,t,u;
  double k,kj,kjp,kjm,logf;
  double sumU,sumL;
  
  h = *tMax / *n;
  t = h;
  u = 0.5*h;
  logf = logrhFct(t,*a,*b);
  Gl[0] = logf;
  G[0] = logf - logkFct(t,u,*a,*b);;
  Gu[0] = logf - logkFct(t,0.0,*a,*b);
  
  for (j = 2; j <= *n; j++) {
    t = j*h;
    sum = logrhFct(t,*a,*b);
    
    for (i = 1; i < j; i++) {
      u=(i-0.5)*h;
      k=logkFct(t,u,*a,*b);
      sum = logspace_sub(sum,k+G[i-1]);
    }
    u=(j-0.5)*h;
    k=logkFct(t,u,*a,*b);
    G[j-1] = sum-k;
  } 
  
  for (i = 1; i < *n; i++) {
    G[i] = logspace_add(G[i-1],G[i]);
  }

  for (i = 1; i < *n; i++) {
    t = (i+1)*h;
    logf = logrhFct(t,*a,*b);
    sumU = logf;
    sumL = logf;
    
    kjm = logkFct(t,0.0,*a,*b);
    kj = logkFct(t,h,*a,*b);
    for (j = 0; j < i; j++) {
      if (j < i-1) {
	u = (j+1)*h;
	kjp = logkFct(t,u+h,*a,*b);
      } else {
	kjp = 0.0;
      }
      sumL = logspace_add(sumL,logspace_sub(kjp,kj) + Gl[j]);
      sumU = logspace_add(sumU,logspace_sub(kj,kjm) + Gu[j]);
      kjm=kj;
      kj=kjp;
    }
    Gl[i] = sumL;
    Gu[i] = sumU - kjm;
  } 
  
  for (i = 0; i < *n; i++) {
    G[i] = exp(G[i]);
    Gl[i] = exp(Gl[i]);
    Gu[i] = exp(Gu[i]);
  }

} 

static double cFct(double t, double alpha, double beta) {
  return(alpha + beta * sqrt(t));
}

static double bFct(double t, double beta) {
  return(0.5*beta/sqrt(t));
}

static double rhFct(double t, double alpha, double beta) {
  double arg1,arg2,arg3,stinv;
  stinv=1/sqrt(t);
  arg1=-cFct(t,alpha,beta)*stinv;
  arg2=-2*bFct(t,beta)*(cFct(t,alpha,beta)-t*bFct(t,beta));
  arg3=(-cFct(t,alpha,beta)+2*t*bFct(t,beta))*stinv;
  return(pnorm(arg1,0.0,1.0,1,0)+exp(arg2)*pnorm(arg3,0.0,1.0,1,0));
}

static double logrhFct(double t, double alpha, double beta) {
  double arg1,arg2,arg3,stinv;
  stinv=1/sqrt(t);
  arg1=-cFct(t,alpha,beta)*stinv;
  arg2=-2*bFct(t,beta)*(cFct(t,alpha,beta)-t*bFct(t,beta));
  arg3=(-cFct(t,alpha,beta)+2*t*bFct(t,beta))*stinv;
  return(logspace_add(pnorm(arg1,0.0,1.0,1,1),arg2+pnorm(arg3,0.0,1.0,1,1)));
}

static double kFct(double t, double u, double alpha, double beta) {
  double arg1,arg2,arg3,sdiffinv;
  sdiffinv=1/sqrt(t-u);
  arg1=(cFct(u,alpha,beta)-cFct(t,alpha,beta))*sdiffinv;
  arg2=-2*bFct(t,beta)*(cFct(t,alpha,beta)-cFct(u,alpha,beta)-(t-u)*bFct(t,beta));
  arg3=(cFct(u,alpha,beta)-cFct(t,alpha,beta)+2*(t-u)*bFct(t,beta))*sdiffinv;
  return(pnorm(arg1,0.0,1.0,1,0)+exp(arg2)*pnorm(arg3,0.0,1.0,1,0));
}

static double logkFct(double t, double u, double alpha, double beta) {
  double arg1,arg2,arg3,sdiffinv;
  sdiffinv=1/sqrt(t-u);
  arg1=(cFct(u,alpha,beta)-cFct(t,alpha,beta))*sdiffinv;
  arg2=-2*bFct(t,beta)*(cFct(t,alpha,beta)-cFct(u,alpha,beta)-(t-u)*bFct(t,beta));
  arg3=(cFct(u,alpha,beta)-cFct(t,alpha,beta)+2*(t-u)*bFct(t,beta))*sdiffinv;
  return(logspace_add(pnorm(arg1,0.0,1.0,1,1),arg2+pnorm(arg3,0.0,1.0,1,1)));
}
