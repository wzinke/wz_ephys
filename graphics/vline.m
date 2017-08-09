function v = vline(X, varargin)
% little hack to easily draw a vertical line without always asking for ylim.
% Just enter the x-coordinate and line properties.
%
% wolf zinke, 17.1.2014

if(nargin<1)
    error('At least the x-coordinate has to be specified!');
end

if(isempty(X))
    return;
end

ax=gca;

if(exist('graph2d.constantline','file'))
    vv = graph2d.constantline(X,'parent',ax);
    vv = changedependvar(vv,'x');
    set(vv,'Color','black','LineWidth',1,'LineStyle','-');
else
    Xp = [X(:), X(:)];
    Yrng = repmat(ylim(),length(X),1);

    vv = plot(Xp, Yrng, 'Color','black','LineWidth',1,'LineStyle','-');
end

set(vv,varargin{:});

if(nargout > 0)
    v=vv;
end



