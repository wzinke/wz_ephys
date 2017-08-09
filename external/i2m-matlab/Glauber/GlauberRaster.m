function raster=GlauberRaster(J1,h,N,M,beta,tau0)
%Model of Glauber Dynamics
if nargin<4
    M=1000000
end
if nargin<5
    beta = 5;
end
if nargin<6
    tau0=1;
end


Vec = rand(N,1);
Vec = 2*ceil(Vec-(1-0.5))-1;

raster=zeros(N,M);

for step=1:M
    VecOld=Vec;
    idx = randperm(N);
    for i=idx
        A=(J1'*Vec + h.*Vec);
        p(i) = (1/(2*tau0))*(1 - Vec(i).*tanh(beta*Vec(i).*A(i)));
        if (rand<p(i))
            Vec(i)=-Vec(i);
        end        
    end
    raster(:,step)=VecOld;
end



