function raster=Generator(h,J,M)
%Monte Carlo evaluation to simulate a surrogate raster
if nargin<3
    M=100000;
end
N=length(h);

Vec = rand(N,1);
Vec = 2*ceil(Vec-0.5)-1;

E = dot(h,Vec)+0.5*dot(Vec,J*Vec);
JV=J*Vec;

for step=1:M
    ind=ceil(N*rand);
    DeltaE = -2*Vec(ind)*h(ind) - 2*Vec(ind)*JV(ind) + 2*J(ind,ind);
    if rand<exp(DeltaE)
        JV=JV-2*Vec(ind)*J(:,ind);
        Vec(ind)=-Vec(ind);
    end
end

 raster = Vec;
 
 