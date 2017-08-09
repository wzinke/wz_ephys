function rast=RasterPopn(r,bin,tmin,tmax)%values must have the same unit (usuallly ms)

a = find(r(1,:)<tmax & r(1,:)>tmin);
s = r(1,a);

s=s - tmin; 

rast = hist(s,round((tmax-tmin)/bin));
