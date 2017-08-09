function wz_figmax(fig)
% maximize a figure to full screen using an undocumented function according to:
% http://undocumentedmatlab.com/blog/minimize-maximize-figure-window/
%
% wolf zinke, 24.1.2014

if(~exist('fig','var') || isempty(fig))
    fig = gcf;
end

jFrame = get(handle(fig),'JavaFrame');
jFrame.setMinimized(true);
