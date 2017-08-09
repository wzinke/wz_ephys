function outcolor = getcolor(varargin)
    

if (~ischar(varargin{1}))
    start(varargin{1},varargin{2})
else
    try
     feval(varargin{:});
    end
end


function start(controlhandle, clusthandles)

global figattrib;
global clustattrib;
global graphattrib;

clustnum = get(controlhandle,'UserData');

if (clustnum > 0)
    initcolor = figattrib.mixcolor(clustnum,:);
    currentcolor = figattrib.mixcolor(clustnum,:);
elseif (clustnum == 0)
    initcolor = clustattrib.cluster0attrib.color;
    currentcolor = clustattrib.cluster0attrib.color;
elseif (clustnum == -1)
    initcolor = graphattrib.backgroundcolor;
    currentcolor = graphattrib.backgroundcolor;    
end

mainfigpos = get(clusthandles.figure1,'Position');

fighandle = figure('MenuBar','none','Color',get(0,'DefaultUicontrolBackgroundColor'),'Position',[round(mainfigpos(1)+(mainfigpos(3)/2)) round(mainfigpos(2)+(mainfigpos(4)/2)) 200 200], ... 
        'Tag','colorChooseFigure','NumberTitle','off','Resize','off','Name',['Cluster ',num2str(clustnum)], ... 
        'CloseRequestFcn','getcolor(''colorChooseFigure_CloseRequestFnc'',gcbo,guidata(gcbo))', ...
        'WindowStyle','modal');

if (clustnum == -1)
    set(fighandle,'Name','Background color');
end
%colaxes = axes('XTick',[],'YTick',[],'Position',[0 0 .2 .2]);
%rectangle('Position',[1 1 25 25],'FaceColor',color);
s(1) = uicontrol('Tag','slider1','Style','slider','Position',[40 50 20 100],'Callback','getcolor(''slider1_Callback'',gcbo,guidata(gcbo))');
s(2) = uicontrol('Tag','slider2','Style','slider','Position',[90 50 20 100],'Callback','getcolor(''slider2_Callback'',gcbo,guidata(gcbo))');
s(3) = uicontrol('Tag','slider3','Style','slider','Position',[140 50 20 100],'Callback','getcolor(''slider3_Callback'',gcbo,guidata(gcbo))');
b(1) = uicontrol('Tag','cancelbutton','Style','pushbutton','Position',[96 5 50 20],'String','Cancel','Callback','getcolor(''cancelbutton_Callback'',gcbo,guidata(gcbo))');
%b(2) = uicontrol('Tag','resetbutton','Style','pushbutton','Position',[96 5 50 20],'String','Default','Callback','getcolor(''resetbutton_Callback'',gcbo,guidata(gcbo))');
b(3) = uicontrol('Tag','okbutton','Style','pushbutton','Position',[147 5 50 20],'String','Ok','Callback','getcolor(''okbutton_Callback'',gcbo,guidata(gcbo))');
t(1) = uicontrol('Tag','text1','Style','text','Position',[30 150 40 20],'String','Red');
t(2) = uicontrol('Tag','text2','Style','text','Position',[80 150 40 20],'String','Green');
t(3) = uicontrol('Tag','text3','Style','text','Position',[130 150 40 20],'String','Blue');
f(1) = uicontrol('Tag','frame1','Style','frame','Position',[20 180 160 10],'BackgroundColor',currentcolor);



set(s(1),'Value',initcolor(1));
set(s(2),'Value',initcolor(2));
set(s(3),'Value',initcolor(3));
handles = guihandles(fighandle);
handles.currentcolor = currentcolor;
handles.initcolor = initcolor;
handles.clustcontrol = clusthandles.clustcontrol;
handles.clustnum = clustnum;
handles.valueok = 0;
handles.fighandle = fighandle;
handles.clusthandles = clusthandles;
guidata(fighandle,handles);
%------------------------------------------------
function slider1_Callback(hObject,handles)
global figattrib;
global clustattrib;
global graphattrib;
clustnum = handles.clustnum;
if (clustnum > 0)
	figattrib.mixcolor(clustnum,1) = get(hObject,'Value');
	set(handles.clustcontrol(clustnum,2),'BackgroundColor',figattrib.mixcolor(handles.clustnum,:));
	set(handles.clustcontrol(clustnum,3),'BackgroundColor',figattrib.mixcolor(handles.clustnum,:));
	set(handles.clustcontrol(clustnum,4),'BackgroundColor',figattrib.mixcolor(handles.clustnum,:));
	set(handles.frame1,'BackgroundColor',figattrib.mixcolor(handles.clustnum,:));
elseif (clustnum == 0)
    clustattrib.cluster0attrib.color(1) = get(hObject,'Value');
    set(handles.clusthandles.cluster0button,'BackgroundColor',clustattrib.cluster0attrib.color)
    set(handles.frame1,'BackgroundColor',clustattrib.cluster0attrib.color);
    if (max(clustattrib.cluster0attrib.color) < .5)
        set(handles.clusthandles.cluster0button,'ForegroundColor',[1 1 1]);
    else
        set(handles.clusthandles.cluster0button,'ForegroundColor',[0 0 0]);
    end
elseif (clustnum == -1)
    graphattrib.backgroundcolor(1) = get(hObject,'Value');
    %set(graphattrib.graphwindow,'Color',graphattrib.backgroundcolor)
    set(handles.frame1,'BackgroundColor',graphattrib.backgroundcolor);
end

%-------------------------------------------------
function slider2_Callback(hObject,handles)
global figattrib;
global clustattrib;
global graphattrib;
clustnum = handles.clustnum;
if (clustnum > 0)
	figattrib.mixcolor(clustnum,2) = get(hObject,'Value');
	set(handles.clustcontrol(clustnum,2),'BackgroundColor',figattrib.mixcolor(handles.clustnum,:));
	set(handles.clustcontrol(clustnum,3),'BackgroundColor',figattrib.mixcolor(handles.clustnum,:));
	set(handles.clustcontrol(clustnum,4),'BackgroundColor',figattrib.mixcolor(handles.clustnum,:));
	set(handles.frame1,'BackgroundColor',figattrib.mixcolor(handles.clustnum,:));
elseif (clustnum == 0)
    clustattrib.cluster0attrib.color(2) = get(hObject,'Value');
    set(handles.clusthandles.cluster0button,'BackgroundColor',clustattrib.cluster0attrib.color)    
    set(handles.frame1,'BackgroundColor',clustattrib.cluster0attrib.color);
    if (max(clustattrib.cluster0attrib.color) < .5)
        set(handles.clusthandles.cluster0button,'ForegroundColor',[1 1 1]);
    else
        set(handles.clusthandles.cluster0button,'ForegroundColor',[0 0 0]);
    end    
elseif (clustnum == -1)    
    graphattrib.backgroundcolor(2) = get(hObject,'Value');
    %set(graphattrib.graphwindow,'Color',graphattrib.backgroundcolor)    
    set(handles.frame1,'BackgroundColor',graphattrib.backgroundcolor);
end
%--------------------------------------------------
function slider3_Callback(hObject,handles)
global figattrib;
global clustattrib;
global graphattrib;
clustnum = handles.clustnum;
if (clustnum > 0)
	figattrib.mixcolor(clustnum,3) = get(hObject,'Value');
	set(handles.clustcontrol(clustnum,2),'BackgroundColor',figattrib.mixcolor(handles.clustnum,:));
	set(handles.clustcontrol(clustnum,3),'BackgroundColor',figattrib.mixcolor(handles.clustnum,:));
	set(handles.clustcontrol(clustnum,4),'BackgroundColor',figattrib.mixcolor(handles.clustnum,:));
	set(handles.frame1,'BackgroundColor',figattrib.mixcolor(handles.clustnum,:));
elseif (clustnum == 0)
    clustattrib.cluster0attrib.color(3) = get(hObject,'Value');
    set(handles.clusthandles.cluster0button,'BackgroundColor',clustattrib.cluster0attrib.color)
    set(handles.frame1,'BackgroundColor',clustattrib.cluster0attrib.color);
    if (max(clustattrib.cluster0attrib.color) < .5)
        set(handles.clusthandles.cluster0button,'ForegroundColor',[1 1 1]);
    else
        set(handles.clusthandles.cluster0button,'ForegroundColor',[0 0 0]);
    end    
elseif (clustnum == -1)    
    graphattrib.backgroundcolor(3) = get(hObject,'Value');
    %set(graphattrib.graphwindow,'Color',graphattrib.backgroundcolor)    
    set(handles.frame1,'BackgroundColor',graphattrib.backgroundcolor);
end
%--------------------------------------------------
function cancelbutton_Callback(hObject,handles)
global figattrib;
global clustattrib;
global graphattrib;
clustnum = handles.clustnum;
if (clustnum > 0)
	figattrib.mixcolor(clustnum,:) = handles.initcolor;
	set(handles.clustcontrol(clustnum,2),'BackgroundColor',figattrib.mixcolor(handles.clustnum,:));
	set(handles.clustcontrol(clustnum,3),'BackgroundColor',figattrib.mixcolor(handles.clustnum,:));
	set(handles.clustcontrol(clustnum,4),'BackgroundColor',figattrib.mixcolor(handles.clustnum,:));
elseif (clustnum == 0)
    clustattrib.cluster0attrib.color = handles.initcolor;
    set(handles.clusthandles.cluster0button,'BackgroundColor',clustattrib.cluster0attrib.color)
    if (max(clustattrib.cluster0attrib.color) < .5)
        set(handles.clusthandles.cluster0button,'ForegroundColor',[1 1 1]);
    else
        set(handles.clusthandles.cluster0button,'ForegroundColor',[0 0 0]);
    end    
elseif (clustnum == -1)    
    graphattrib.backgroundcolor = handles.initcolor;
    %set(graphattrib.graphwindow,'Color',graphattrib.backgroundcolor)
end
matclust('plotgraph',handles.clusthandles);
delete(handles.fighandle);
%--------------------------------------------------
function resetbutton_Callback(hObject,handles)
global figattrib;
global clustattrib;
global graphattrib;
load matclust_defaults;
clustnum = handles.clustnum;
if (clustnum > 0)
	figattrib.mixcolor(clustnum,:) = ClusterColors(clustnum,:);
	set(handles.slider1,'Value',ClusterColors(clustnum,1));
	set(handles.slider2,'Value',ClusterColors(clustnum,2));
	set(handles.slider3,'Value',ClusterColors(clustnum,3));
	set(handles.clustcontrol(clustnum,2),'BackgroundColor',figattrib.mixcolor(handles.clustnum,:));
	set(handles.clustcontrol(clustnum,3),'BackgroundColor',figattrib.mixcolor(handles.clustnum,:));
	set(handles.clustcontrol(clustnum,4),'BackgroundColor',figattrib.mixcolor(handles.clustnum,:));
	set(handles.frame1,'BackgroundColor',figattrib.mixcolor(handles.clustnum,:));
elseif (clustnum == 0)
    clustattrib.cluster0attrib.color = Cluster0Color;
    set(handles.slider1,'Value',clustattrib.cluster0attrib.color(1));
	set(handles.slider2,'Value',clustattrib.cluster0attrib.color(2));
	set(handles.slider3,'Value',clustattrib.cluster0attrib.color(3));
    set(handles.clusthandles.cluster0button,'BackgroundColor',clustattrib.cluster0attrib.color);
    set(handles.frame1,'BackgroundColor',clustattrib.cluster0attrib.color);
    if (max(clustattrib.cluster0attrib.color) < .5)
        set(handles.clusthandles.cluster0button,'ForegroundColor',[1 1 1]);
    else
        set(handles.clusthandles.cluster0button,'ForegroundColor',[0 0 0]);
    end    
elseif (clustnum == -1)    
    graphattrib.backgroundcolor = GraphBackgroundColor;
    set(handles.slider1,'Value',graphattrib.backgroundcolor(1));
	set(handles.slider2,'Value',graphattrib.backgroundcolor(2));
	set(handles.slider3,'Value',graphattrib.backgroundcolor(3));
    %set(graphattrib.graphwindow,'Color',graphattrib.backgroundcolor);
    set(handles.frame1,'BackgroundColor',graphattrib.backgroundcolor);
end
%--------------------------------------------------
function okbutton_Callback(hObject,handles)
global figattrib;
global graphattrib;
global clustattrib;

clustnum = handles.clustnum;
handles.valueok = 1;
guidata(handles.colorChooseFigure,handles);
currdir = pwd;
cd(figattrib.foldername);
load matclust_defaults;




%clustattrib.cluster0attrib.color = Cluster0Color; %stores the color of cluster 0
%graphattrib.backgroundcolor = GraphBackgroundColor; %stores the color of the graphs background
%figattrib.mixcolor = ClusterColors; %stores the colors of the clusters
if (clustnum == -1)
    set(graphattrib.graphwindow,'Color',graphattrib.backgroundcolor);
    matclust_defaults.GraphBackgroundColor = graphattrib.backgroundcolor;
elseif (clustnum == 0)
    matclust_defaults.Cluster0Color = clustattrib.cluster0attrib.color;
elseif (clustnum > 0)
    matclust_defaults.ClusterColors = figattrib.mixcolor;
    set(figattrib.clustcontrol(clustnum,1),'BackgroundColor',figattrib.mixcolor(clustnum,:));
end



save matclust_defaults matclust_defaults;
cd(currdir);

matclust('plotgraph',handles.clusthandles);
delete(handles.fighandle);
%-------------------------------------------------
function colorChooseFigure_CloseRequestFnc(hObject,handles)
global figattrib;
global clustattrib;
global graphattrib;

response = questdlg('Keep color selection?','MATCLUST');
if (strcmp(response,'No'))
    clustnum = handles.clustnum;
    if (clustnum > 0)
        figattrib.mixcolor(clustnum,:) = handles.initcolor;
        set(handles.clustcontrol(clustnum,2),'BackgroundColor',figattrib.mixcolor(handles.clustnum,:));
        set(handles.clustcontrol(clustnum,3),'BackgroundColor',figattrib.mixcolor(handles.clustnum,:));
        set(handles.clustcontrol(clustnum,4),'BackgroundColor',figattrib.mixcolor(handles.clustnum,:));
    elseif (clustnum == 0)
        clustattrib.cluster0attrib.color = handles.initcolor;
        set(handles.clusthandles.cluster0button,'BackgroundColor',clustattrib.cluster0attrib.color)
        if (max(clustattrib.cluster0attrib.color) < .5)
            set(handles.clusthandles.cluster0button,'ForegroundColor',[1 1 1]);
        else
            set(handles.clusthandles.cluster0button,'ForegroundColor',[0 0 0]);
        end
    elseif (clustnum == -1)    
        graphattrib.backgroundcolor = handles.initcolor;
        %set(graphattrib.graphwindow,'BackgroundColor',graphattrib.backgroundcolor)
    end
    delete(hObject)
elseif (strcmp(response,'Yes'))
    if (clustnum == -1)
        set(graphattrib.graphwindow,'Color',graphattrib.backgroundcolor);
    end
    matclust('plotgraph',handles.clusthandles);
    delete(hObject)      
end

