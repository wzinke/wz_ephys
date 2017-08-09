function Maxima = maxima(x,boolExtremes)
% Find location of local maxima excluding extremes

maxflags=find(diff(sign(diff([-inf;x(:);-inf])))<0);

% Getting rid of extremes
if(boolExtremes)
    Maxima=maxflags(maxflags~=1 & maxflags~=length(x));
else
    Maxima=maxflags;
end
