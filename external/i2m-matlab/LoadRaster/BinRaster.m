function rast=BinRaster(r,bin,tmin,tmax)%values must have the same unit (usually ms)

a = find(r(1,:)<tmax & r(1,:)>tmin);
s = r(:,a);

nbneur=max(s(2,:))

s(1,:)=floor((s(1,:) - min(s(1,:)) )/bin) +1; 

len=max(s(1,:));

rast=zeros(nbneur,len);

rast(:,:)=-1;

list=nbneur*(s(1,:)-1)+1;
list=list + s(2,:) - 1;

rast(list)=1;

