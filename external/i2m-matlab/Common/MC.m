function [m,C]=MC(h,J,N,M)
%Monte Carlo evaluation
if nargin<3
    M=100000
end

Vec = rand(N,1);
Vec = 2*ceil(Vec-0.5)-1;

E=dot(h,Vec)+0.5*dot(Vec,J*Vec);
JV=J*Vec;

sum_m=zeros(N,1);
sum_C=zeros(N,N);

for step=1:M
    ind=ceil(N*rand);
    DeltaE = -2*Vec(ind)*h(ind) - 2*Vec(ind)*JV(ind) + 2*J(ind,ind);
    if rand<exp(DeltaE)
        JV=JV-2*Vec(ind)*J(:,ind);
        Vec(ind)=-Vec(ind);
    end
    sum_m=sum_m+Vec;
    sum_C=sum_C+Vec*transpose(Vec);
end

m=sum_m/M;
C=sum_C/M;
