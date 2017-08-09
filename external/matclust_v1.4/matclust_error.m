function matclust_error(varargin)

if (nargin == 0)
    start;
else
    feval(varargin{:});   
end

%--------------------------------------------------------------------
function start
%called to create the figure


w = lasterr;
returnchar = findstr(w,char(10));
if isempty(returnchar)
    returnchar = 0;
end
w = w(returnchar+1:end);
%create the figure
fighandle = errordlg(w,'Error');


%create the widgets

 

b(1) = uicontrol(fighandle,'Tag','infobutton','Style','pushbutton','Units','pixels','Position',[3 3 15 15],'String','i','Callback','matclust_error(''infobutton_Callback'',guidata(gcbo))');

%attach important variables to the figure handle
fighandles = guihandles(fighandle);

guidata(fighandle,fighandles);
%------------------------------------------------------------
function infobutton_Callback(handles)

w = lasterr;
returnchar = findstr(w,char(10));
if isempty(returnchar)
    returnchar = length(w)+1;
end
w = w(1:returnchar-1);
msgbox(w);