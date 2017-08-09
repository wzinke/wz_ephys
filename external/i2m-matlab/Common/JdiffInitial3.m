function diffres=JdiffInitial3(J,m,A)
%Intermediate function to minimize to compute the third order approximation
%of J

diffres=0;

N=size(J,1);

for i=1:N
    for j=1:N
        if j~=i
            difftemp=A(i,j)+J(i,j)+2*(J(i,j)^2)*m(i)*m(j)...
                +(2/3)*(J(i,j)^3)*(1-3*m(i)^3)*(1-3*m(j)^3);
            for k=1:N
                if (k~=i)&&(k~=j)
                    difftemp=difftemp+4*J(i,j)*J(i,k)*J(j,k)*m(i)*m(j)*(1-m(k)^2);
                end
            end
            diffres=diffres+difftemp^2;
        end
    end
end

diffres = sqrt(diffres/(N*(N-1)));
