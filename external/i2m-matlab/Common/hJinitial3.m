function [h,J] = hJinitial3(m,C,N)
%Approximation taken from Tanaka et al
A=inv(C-m*transpose(m));

[h2,J2] = hJinitial2(m,C,N);

[J,diffres]=fminsearch(@(Jvar) JdiffInitial3(Jvar,m,A),J2,optimset('Display','notify'));%,'TolFun',1e-9,'TolX',1e-9));

for i=1:N
    J(i,i)=1/(1-m(i)*m(i))-A(i,i);
end

h=zeros(N,1);

for i=1:N
    h(i)=atanh(m(i));
    for j=1:N
        if (j~=i)
            h(i)= h(i) -J(i,j)*m(j)*(1+(2/3)*(J(i,j)^2)*(1+3*m(i)^2)*(1-m(j)^2)) ;
        else
            h(i) = h(i) -J(i,j)*m(j);
        end
    end
end

