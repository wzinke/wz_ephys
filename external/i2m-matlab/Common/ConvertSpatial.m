function [mtot,Ctot] = ConvertSpatial(m,C,C1)
%Create the vector and the matrix for the 2*N spin vector
N=length(m);
mtot(1:N,1)=m;
mtot((N+1):2*N,1)=m;

Ctot(1:N , 1:N) = C;
Ctot(1:N , (N+1):2*N) = C1;
Ctot((N+1):2*N , 1:N) = C1';
Ctot((N+1):2*N , (N+1):2*N) = C;

