function otherfilteradd(varargin)

if (~ischar(varargin{1}))
    start(varargin{1},varargin{2})
else
     feval(varargin{:});
     
end


function start(hObject,handles)

mainfighandles = handles;

mainfigpos = get(handles.figure1,'Position');

fighandle = figure('MenuBar','none','Color',get(0,'DefaultUicontrolBackgroundColor'),'Position',[round(mainfigpos(1)+(mainfigpos(3)/2)) round(mainfigpos(2)+(mainfigpos(4)/2)) 230 200], ... 
        'Tag','addfilterFigure','NumberTitle','off','Resize','off','Name','Add filter', ... 
        'CloseRequestFcn','otherfilteradd(''addfilterFigure_CloseRequestFnc'',gcbo,guidata(gcbo))', ...
        'WindowStyle','modal');

e(1) = uicontrol('Tag','edit1','Style','edit','Position',[70 152 115 20],'BackgroundColor',[1 1 1]);
e(2) = uicontrol('Tag','edit2','Style','edit','Position',[70 102 115 20],'BackgroundColor',[1 1 1]);
%e(3) = uicontrol('Tag','edit3','Style','edit','Position',[70 52 115 20],'BackgroundColor',[1 1 1]);
b(1) = uicontrol('Tag','cancelbutton','Style','pushbutton','Position',[126 5 50 20],'String','Cancel','Callback','otherfilteradd(''cancelbutton_Callback'',gcbo,guidata(gcbo))');
b(2) = uicontrol('Tag','okbutton','Style','pushbutton','Position',[177 5 50 20],'String','Ok','Callback','otherfilteradd(''okbutton_Callback'',gcbo,guidata(gcbo))');
b(3) = uicontrol('Tag','lookbutton','Style','pushbutton','Position',[187 102 25 20],'String','/','Callback','otherfilteradd(''lookbutton_Callback'',gcbo,guidata(gcbo))');
t(1) = uicontrol('Tag','text1','Style','text','Position',[15 150 50 20],'String','Name:','HorizontalAlignment','left');
t(2) = uicontrol('Tag','text2','Style','text','Position',[15 100 50 20],'String','Program:','HorizontalAlignment','left');
%t(3) = uicontrol('Tag','text3','Style','text','Position',[15 50 50 20],'String','End time:','HorizontalAlignment','left');
t(4) = uicontrol('Tag','text4','Style','text','Position',[5 180 200 20],'String','','HorizontalAlignment','left','ForegroundColor',[1 0 0],'BackgroundColor',get(0,'DefaultUicontrolBackgroundColor'));
fighandles = guihandles(fighandle);
fighandles.fighandle = fighandle;
fighandles.mainfighandles = mainfighandles;
fighandles.pathname = pwd;
guidata(fighandle,fighandles);
%------------------------------------------------------------
function addfilterFigure_CloseRequestFnc(hObject,fighandles)

delete(hObject);
%------------------------------------------------------------
function cancelbutton_Callback(hObject,fighandles)

delete(fighandles.addfilterFigure);
%------------------------------------------------------------
function lookbutton_Callback(hObject,fighandles)
global figattrib

currdir = pwd;
cd(figattrib.foldername);
cd('Filters');

if (ispc)
    [filename, pathname] = uigetfile('*.m','Find filter program');
else
    [filename, pathname] = filebrowse('open','filter','*.m','title','Find filter program');
end
if ischar(filename)
    set(fighandles.edit2,'String',[filename]);
    fighandles.pathname = pathname;
end
guidata(fighandles.fighandle,fighandles);

cd(currdir);
%---------------------------------------------------------

function okbutton_Callback(hObject,fighandles)

global clustdata;


origpath = pwd;
eval(['cd ',fighandles.pathname]);


names = get(fighandles.mainfighandles.otherfilterList,'String');
currval = get(fighandles.mainfighandles.otherfilterList,'Value');

set(fighandles.text4,'String','Calculating...');
drawnow
%try
    pname = get(fighandles.edit2,'String');
    pointfind = strfind(pname,'.');
    if ~isempty(pointfind)
        pname = pname(1:pointfind-1);
    end
    %pwd
    outfilt = feval(pname);
% catch
%     set(fighandles.text4,'String',[pname,' is not a valid filter.  yo yo']);
%     cd(pwd)
%     return        
% end

if ((~isequal(size(outfilt),[size(clustdata.params,1),1]))|(max(outfilt)>1))
    set(fighandles.text4,'String',[pname,' is not a valid filter.']);
    cd(origpath)
    return    
end
    
eval(['cd ',origpath]);
newname = [get(fighandles.edit1,'String')];
newname = [num2str(currval),'  ',newname];

names{currval} = newname;
%clustdata.timefiltermemmap(find(clustdata.timefiltermemmap>currval)) = clustdata.timefiltermemmap(find(clustdata.timefiltermemmap>currval))+1;
    
outfilt = logical(outfilt);
set(fighandles.mainfighandles.otherfilterList,'String',names);
clustdata.otherfilternames = names;
%firstzero = min(find(clustdata.otherfiltermemmap==0));
%clustdata.otherfiltermemmap(firstzero) = currval;
%clustdata.otherfilters = fastbitset(clustdata.otherfilters,firstzero,outfilt);
clustdata.otherfiltermemmap(currval) = currval;
clustdata.otherfilters = fastbitset(clustdata.otherfilters,currval,outfilt);


%set(fighandles.mainfighandles.timefilterList,'Value',[]);
matclust('clearotherfilter');
clustdata.filtermemmap = [clustdata.timefiltermemmap; clustdata.otherfiltermemmap];

matclust('addnewstate','add filter',fighandles.mainfighandles); 
delete(fighandles.addfilterFigure);