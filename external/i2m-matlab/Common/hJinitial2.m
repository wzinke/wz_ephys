function [h,J] = hJinitial2(m,C,N)
%Approximation taken from Tanaka et al
A=inv(C-m*transpose(m));

h=zeros(N,1);

for i=1:N
    for j=i+1:N
        J(i,j)=real(-1+sqrt(1-8*A(i,j)*m(i)*m(j)))/(4*m(i)*m(j));
        J(j,i)=J(i,j);
    end
    J(i,i)=1/(1-m(i)*m(i))-A(i,i);
end

for i=1:N
    h(i)=atanh(m(i));
    for j=1:N
        h(i)= h(i) -J(i,j)*m(j);
    end
end
