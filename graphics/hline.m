function h = hline(Y, varargin)
% little hack to easily draw a horizontal line without always asking for xlim.
% Just enter the y-coordinate and line properties.
%
% wolf zinke, 17.1.2014

if(nargin<1)
    error('At least the y-coordinate has to be specified!');
end

if(isempty(Y))
    return;
end

ax=gca;

if(exist('graph2d.constantline','file'))
    hh = graph2d.constantline(Y,'parent',ax);
    hh = changedependvar(hh,'y');
    set(hh,'Color','black','LineWidth',1,'LineStyle','-');
else
    Yp = [Y(:), Y(:)];
    Xrng = repmat(xlim(),1,length(Y));
    hh = plot(Xrng, Yp, 'Color','black','LineWidth',1,'LineStyle','-');
end
set(hh,varargin{:});

if(nargout > 0)
    h=hh;
end


