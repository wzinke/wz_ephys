function [h,J,J1,hr,Jr] = DeConvertSpatial(htot,Jtot,hstat,Jstat)
N=length(htot)/2;
h=htot(N+1:2*N,1);
hr=h(1:N,1)-hstat;
Jr=Jtot(1:N , 1:N)-Jstat;

J=Jtot(N+1:2*N , N+1:2*N);
J1=0.5*( Jtot(1:N , (N+1):2*N) + Jtot( (N+1):2*N , 1:N)' );

