function fig = nice_plot(varargin)

nm = 'my nice figure';
figbck = [1,1,1];
fig_pos = [0 0.05 1 0.85];

if(nargin > 0)
    for(i=1:nargin)
        if(ischar(varargin{i}) == 1)
            nm = varargin{i};
        elseif(isnumeric(varargin{i}) == 1)
            if(length(varargin{i}) == 1)
                fig = Varargin{i};
            elseif(length(varargin{i}) == 4)                
                figbck = varargin{i};
            elseif(length(varargin{i}) == 4)
                fig_pos = varargin{i};
            end
        end
    end
end

if(exist('fig','var') == 0 || isempty(fig) == 1)
    fig = gcf;
end

set(fig, 'Name', nm);
set(fig, 'Color',figbck);

set(gca,'fontsize',14);

set(fig, 'Units', 'Normalized');
set(fig, 'Renderer', 'Painters ');

set(fig, 'NumberTitle', 'off');

set(gca,'TickDir','out');

axis tight

% set(fig,'PaperOrientation','landscape');
% % set(fig,'PaperOrientation','portrait');
% set(fig,'PaperUnits','normalized');
% set(fig,'PaperType','A4');
% % set(fig,'PaperSize',[297, 210]);
% set(fig,'PaperPosition',[0,0,1,1]);
