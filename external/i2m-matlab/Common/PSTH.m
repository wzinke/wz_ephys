function [rast,x]=PSTH(r,bin,tmin,tmax)%values must have the same unit (usuallly ms)

if nargin<3st,
    tmin = min(r(1,:))-1;
    tmax = max(r(1,:))+1;
end

a = find(r(1,:)<tmax & r(1,:)>tmin);
s = r(:,a);

nbneur=max(s(2,:))

s(1,:)=floor((s(1,:) - min(s(1,:)) )/bin) +1; 

len=max(s(1,:));

rast=zeros(nbneur,len);

list=nbneur*(s(1,:)-1)+1;
list=list + s(2,:) - 1;

for i=1:length(list)
    rast(list(i))=rast(list(i)) + 1;
end

x=(1:size(rast,2))*bin;

