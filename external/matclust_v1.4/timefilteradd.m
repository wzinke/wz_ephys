function timefilteradd(varargin)

if (~ischar(varargin{1}))
    start(varargin{1},varargin{2})
else
     feval(varargin{:});
     
end


function start(hObject,handles)

mainfighandles = handles;

mainfigpos = get(handles.figure1,'Position');

fighandle = figure('MenuBar','none','Color',get(0,'DefaultUicontrolBackgroundColor'),'Position',[round(mainfigpos(1)+(mainfigpos(3)/2)) round(mainfigpos(2)+(mainfigpos(4)/2)) 200 200], ... 
        'Tag','addtimeFigure','NumberTitle','off','Resize','off','Name','Add time filter', ... 
        'CloseRequestFcn','timefilteradd(''addtimeFigure_CloseRequestFnc'',gcbo,guidata(gcbo))', ...
        'WindowStyle','modal');

e(1) = uicontrol('Tag','edit1','Style','edit','Position',[70 152 115 20],'BackgroundColor',[1 1 1]);
e(2) = uicontrol('Tag','edit2','Style','edit','Position',[70 102 115 20],'BackgroundColor',[1 1 1]);
e(3) = uicontrol('Tag','edit3','Style','edit','Position',[70 52 115 20],'BackgroundColor',[1 1 1]);
b(1) = uicontrol('Tag','cancelbutton','Style','pushbutton','Position',[96 5 50 20],'String','Cancel','Callback','timefilteradd(''cancelbutton_Callback'',gcbo,guidata(gcbo))');
b(2) = uicontrol('Tag','okbutton','Style','pushbutton','Position',[147 5 50 20],'String','Ok','Callback','timefilteradd(''okbutton_Callback'',gcbo,guidata(gcbo))');
t(1) = uicontrol('Tag','text1','Style','text','FontUnits','pixels','FontSize',10,'Position',[15 150 50 20],'String','Name:','HorizontalAlignment','left');
t(2) = uicontrol('Tag','text2','Style','text','FontUnits','pixels','FontSize',10,'Position',[15 100 50 20],'String','Start time:','HorizontalAlignment','left');
t(3) = uicontrol('Tag','text3','Style','text','FontUnits','pixels','FontSize',10,'Position',[15 50 50 20],'String','End time:','HorizontalAlignment','left');
t(4) = uicontrol('Tag','text4','Style','text','FontUnits','pixels','FontSize',10,'Position',[5 180 200 20],'String','','HorizontalAlignment','left','ForegroundColor',[1 0 0],'BackgroundColor',get(0,'DefaultUicontrolBackgroundColor'));
fighandles = guihandles(fighandle);
fighandles.mainfighandles = mainfighandles;
guidata(fighandle,fighandles);
%------------------------------------------------------------
function addtimeFigure_CloseRequestFnc(hObject,fighandles)

delete(hObject);
%------------------------------------------------------------
function cancelbutton_Callback(hObject,fighandles)

delete(fighandles.addtimeFigure);
%------------------------------------------------------------
function okbutton_Callback(hObject,fighandles)

global clustdata;

starttime{1} = get(fighandles.edit2,'String');
endtime{1} = get(fighandles.edit3,'String');
timerange = clustdata.datarange(:,1);
names = get(fighandles.mainfighandles.timefilterList,'String');
currval = get(fighandles.mainfighandles.timefilterList,'Value');
try
    t1 = timetrans(starttime,clustdata.UnitsPerSec,2);
    t2 = timetrans(endtime,clustdata.UnitsPerSec,2);
catch
    set(fighandles.text4,'String','Time format- hours:minutes:seconds');
    set(fighandles.edit2,'String','');
    set(fighandles.edit3,'String','');
    return
end
if ((t1>=timerange(1))&(t1<=timerange(2))&(t2>=timerange(1))&(t2<=timerange(2))&(t2>t1))
    newname = [get(fighandles.edit1,'String'),' ',get(fighandles.edit2,'String'),'-',get(fighandles.edit3,'String')];
    newname = [num2str(currval),'  ',newname];
    
    names{currval} = newname;
    %clustdata.timefiltermemmap(find(clustdata.timefiltermemmap>currval)) = clustdata.timefiltermemmap(find(clustdata.timefiltermemmap>currval))+1;
        
    clustdata.timefilterranges(currval,1:2) = [t1 t2];
    set(fighandles.mainfighandles.timefilterList,'String',names);
    clustdata.timefilternames = names;
    passfilter = ((clustdata.params(:,1) >= t1)&(clustdata.params(:,1) <= t2));
    %firstzero = min(find(clustdata.timefiltermemmap==0));
    %clustdata.timefiltermemmap(firstzero) = currval;
    %clustdata.timefilters = fastbitset(clustdata.timefilters,firstzero,passfilter);
    clustdata.timefiltermemmap(currval) = currval;
    clustdata.timefilters = fastbitset(clustdata.timefilters,currval,passfilter);
    %set(fighandles.mainfighandles.timefilterList,'Value',[]);
    matclust('cleartimefilter');
    clustdata.filtermemmap = [clustdata.timefiltermemmap; clustdata.otherfiltermemmap];
else
    trange = [clustdata.params(1,1);clustdata.params(end,1)];
    trange = timetrans(trange,clustdata.UnitsPerSec,1);
    set(fighandles.text4,'String',['Time range- ', trange{1},' to ',trange{2}]);
    set(fighandles.edit2,'String','');
    set(fighandles.edit3,'String','');
    return
end  
matclust('addnewstate','add time filter',fighandles.mainfighandles); 
delete(fighandles.addtimeFigure);