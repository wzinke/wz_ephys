function p = mypcolor(varargin)
% mypcolor  simple wrapper for pcolor
% 
% Unlike the built-in pcolor command, this function does not cut off the last
% last row and column of the input matrix, neither it works wrong on non-regular
% axis.
% 
% Usage: p = mypcolor([x,y],z,[logx],[logy]);
%
%INPUTS:
%
%    x: an array of x values.  can also specify a min and max x value.
%       these values correspond to columns of z. [IF THIS ARGUMENT IS USED,
%       MUST ALSO SPECIFY Y VALUES.]
% 
%    y: an array of y values.  can also specify a min and max y value.
%       these values correspond to rows of z.  [IF THIS ARGUMENT IS USED,
%       MUST ALSO SPECIFY X VALUES.]
% 
%    z: a 2d matrix of values.  this matrix determines the color at each
%       point.
% 
% logx: if this optional argument is set to true, the x-axis will plotted
%       in log scale (similar to semilogx).
%
% logy: if this optional argument is set to true, the y-axis will plotted
%       in log scale (similar to semilogy).
%
%OUTPUTS:
%
%    p: a handle to the resulting pcolor image.
%
% EXAMPLE:
%
%   m = membrane;
%   p = mypcolor(m);
%
% SEE ALSO: PCOLOR, IMAGE, IMAGESC, SEMILOGX, SEMILOGY, LOGLOG, PADARRAY
%
%   AUTHOR: Salva Ardid
%   CONTACT: sardid@yorku.ca
%   BASED ON sanepcolor by Jeremy R. Manning


% CHANGELOG
% 30/06/2012    changes in sanepcolor to avoid problems on non-regular axis

%parse arguments
if length(varargin) == 1 %just a z
    z = varargin{1};
    x = [1 size(z,2)];
    y = [1 size(z,1)];
    [logx,logy] = deal(false);    
elseif (length(varargin) >= 4) %x, y, z, and possibly logx and/or logy
    x = varargin{1};
    y = varargin{2};
    z = varargin{3};
    logx = varargin{4};
    if length(varargin) >= 5, logy = varargin{5}; else logy = false; end
elseif length(varargin) == 2 %z and logx
    z = varargin{1};
    logx = varargin{2};
    logy = false;
    x = [1 size(z,2)];
    y = [1 size(z,1)];
else %length(varargin) == 3
    if isempty(varargin)
        fprintf('\nUsage: p = sanePColor([x,y],z,[logx],[logy]);\n');
        fprintf('Type ''help %s'' for more info.\n\n',mfilename);
        p = [];
        return;
    end
    %posibility 1: x, y, z
    if length(varargin{1}) >= 1 && length(varargin{2}) >= 1
        x = varargin{1};
        y = varargin{2};
        z = varargin{3};
        [logx,logy] = deal(false);
    %posibility 2: z, logx, and logy
    else
        z = varargin{1};
        x = [1 size(z,2)];
        y = [1 size(z,1)];
        logx = varargin{2};
        logy = varargin{3};
    end
end

z = padarray(z,[1 1],'replicate','post');
% following changes are to make it work with non-regular axis
x=x(:);
y=y(:);

if logx
    %newx = logspace(log10(min(x)),log10(max(x)),size(z,2));
    if(length(x)>1)
        newx = [log10(x);log10(x(end)+(x(end)-x(end-1)))];
    else
        newx = [log10(x);log10(x(end)+1)];
    end
else
    %newx = linspace(min(x),max(x),size(z,2));
    if(length(x)>1)
        newx = [x;x(end)+(x(end)-x(end-1))];
    else
        newx = [x;x(end)+1];
    end
end

if logy
    %newy = logspace(log10(min(y)),log10(max(y)),size(z,1));
    if(length(y)>1)
        newy = [log10(y);log10(y(end)+(y(end)-y(end-1)))];
    else
        newy = [log10(y);log10(y(end)+1)];
    end
else
    %newy = linspace(min(y),max(y),size(z,1));
    if(length(y)>1)
        newy = [y;y(end)+(y(end)-y(end-1))];
    else
        newy = [y;y(end)+1];
    end
end

p = pcolor(newx,newy,z);
set(gca,'XTickMode','Manual','YTickMode','Manual');
if exist('logx','var') && logx
    set(gca,'XScale','log');
    if(length(x)>1)
        set(gca,'XTick',logspace(log10(min(x)),log10(max(x)),4));
    end
else
    if(length(x)>1)
        set(gca,'XTick',linspace(min(x),max(x),4));
    end
end
if exist('logy','var') && logy
    set(gca,'YScale','log');
    if(length(y)>1)
        set(gca,'YTick',logspace(log10(min(y)),log10(max(y)),8));
    end
else
    if(length(y)>1)
        set(gca,'YTick',linspace(min(y),max(y),8));
    end
end

shading flat;