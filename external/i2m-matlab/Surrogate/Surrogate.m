function raster=Surrogate(h,J,J1,datalen,M)
if nargin<5
    M=1000
end

N=length(h);

Vec = rand(N,1);
Vec = 2*ceil(Vec-0.5)-1;

raster=zeros(N,datalen);

for step=1:datalen
    raster(:,step)=Vec;
    Vec=Generator((h+J1'*Vec),J,M);
end



