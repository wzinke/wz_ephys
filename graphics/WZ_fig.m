function fig = WZ_fig(varargin)

fig =  figure;

nm = 'WZ standard figure';
figbck = [1,1,1];
fig_pos = [0 0.05 1 0.85];

if(nargin > 0)
    for(i=1:nargin)
        if(ischar(varargin{i}) == 1)
            nm = varargin{i};
        elseif(isnumeric(varargin{i}) == 1)
            if(length(varargin{i}) == 3)
                figbck = varargin{i};
            elseif(length(varargin{i}) == 4)
                fig_pos = varargin{i};
            end
        end
    end
end

set(fig, 'Name', nm);
set(fig, 'Color',figbck);

set(fig, 'NumberTitle', 'off');
set(fig,'PaperOrientation','landscape');
% set(fig,'PaperOrientation','portrait');
set(fig,'PaperUnits','normalized');
set(fig,'PaperType','A4');
% set(fig,'PaperSize',[297, 210]);
set(fig,'PaperPosition',[0,0,1,1]);

set(fig, 'Renderer', 'Painters ');
% set(fig, 'Renderer', 'OpenGL');
set(fig, 'BackingStore', 'on');
set(fig, 'DoubleBuffer', 'on');


warning off MATLAB:HandleGraphics:RenamedProperty:XTickLabels
warning off MATLAB:HandleGraphics:RenamedProperty:YTickLabels
clf;
