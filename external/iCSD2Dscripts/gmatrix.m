function G=gmatrix (n,method)

% GMATRIX Calculates the matrix G.
%   Helper script for two-dimensional ICSD.

% G=gmatrix (n,method)
% Calculates the matrix G 
% which defines the system of equations relating
% the second derivative to the values of the interpolating spline.
% n - number of nodes,
% method = 'n' (natural splines) or 'm' (not-a-knot splines)

%   Copyright 2006-2010 Szymon Leski
%   s.leski@nencki.gov.pl

if method=='n' 
    natural_splines = 1; 
end;
if method=='m' 
    natural_splines = 0; 
end;

a1 = zeros(n-1,1);
a1(1:n-2)=1/6;
a2 = ones(n,1);
a2(2:n-1)=2/3;
a3 = zeros(n-1,1);
a3(2:n-1)=1/6;
A=diag(a1,-1)+diag(a2,0)+diag(a3,1);
if natural_splines==0
    A(n,n-2) = 1;
    A(n,n-1) = -2;
    A(1,2) = -2;
    A(1,3)=1;
end;

b1 = ones(n-1,1);
b1(n-1)=0;
b2 = -2*ones(n,1);
b2(1)=0;
b2(n)=0;
b3 = ones(n-1,1);
b3(1)=0;
B=diag(b1,-1)+diag(b2,0)+diag(b3,1);

G=inv(A)*B;
