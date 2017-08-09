function [h,x,y] = ellipse(x,y,r, varargin)

if(length(r) == 1)
    hr = r;
    vr = r;
else
    hr = r(1);
    vr = r(2);    
end

ang=-pi:0.01:pi;
x=x+hr*cos(ang);
y=y+vr*sin(ang);

hh = plot(x,y, 'color','k','LineStyle','-');    

set(hh,varargin{:});

if(nargout > 0)
    h=hh;
end

