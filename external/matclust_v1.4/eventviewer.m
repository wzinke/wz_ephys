function eventviewer(varargin)

if (~ischar(varargin{1}))
    start(varargin{1},varargin{2})
else
     feval(varargin{:});
     
end

%--------------------------------------------------------------------
function start(clustnum,handles)
%called to create the figure

global clustdata;
global clustattrib;
global figattrib;

mainfighandles = handles;

mainfigpos = get(handles.figure1,'Position');

%create the figure
fighandle = figure('MenuBar','none','Color',get(0,'DefaultUicontrolBackgroundColor'),'Position',[round(mainfigpos(1)+(mainfigpos(3)/2)) round(mainfigpos(2)) 500 500], ... 
        'Tag','eventviewer','NumberTitle','off','Resize','off','Name',['Spike viewer: Cluster ',num2str(clustnum)], ... 
        'CloseRequestFcn','eventviewer(''eventviewer_CloseRequestFnc'',gcbo,guidata(gcbo))', ...
        'WindowStyle','normal');

set(fighandle,'KeyPressFcn','eventviewer(''eventviewer_KeyPressFcn'',gcbo,guidata(gcbo))');
set(fighandle,'WindowButtonDownFcn','eventviewer(''eventviewer_WindowButtonDownFcn'',gcbo,guidata(gcbo))');
set(fighandle,'WindowButtonUpFcn','eventviewer(''eventviewer_WindowButtonUpFcn'',gcbo,guidata(gcbo))');
set(fighandle,'WindowButtonMotionFcn','eventviewer(''eventviewer_WindowButtonMotionFcn'',gcbo,guidata(gcbo))');

b = get(fighandle,'Position');   
clustercolor = figattrib.mixcolor(clustnum,:);
%create the widgets

f(1) = uicontrol('Style','frame','Tag','frame1','Units','normalized','BackgroundColor',get(0,'DefaultUicontrolBackgroundColor'),'ForegroundColor', [0 0 0],'Position',[.01 .52 .98 .15]);
f(2) = uicontrol('Style','frame','Tag','frame2','Units','normalized','BackgroundColor',get(0,'DefaultUicontrolBackgroundColor'),'ForegroundColor', [0 0 0],'Position',[.52 .1 .47 .41]);
f(3) = uicontrol('Style','frame','Tag','frame3','Units','normalized','BackgroundColor',get(0,'DefaultUicontrolBackgroundColor'),'ForegroundColor', [0 0 0],'Position',[.52 .01 .47 .08]);
f(4) = uicontrol('Style','frame','Tag','frame4','Units','normalized','BackgroundColor',get(0,'DefaultUicontrolBackgroundColor'),'ForegroundColor', [0 0 0],'Position',[.01 .52 .43 .15]);

q = [.7 .7 .7];  %axes background color
a(1) = axes('Tag','axes1','Position',[.01 .7 .23 .25],'XTick',[],'YTick',[],'NextPlot','replacechildren','Color',q(1,:));
a(2) = axes('Tag','axes2','Position',[.26 .7 .23 .25],'XTick',[],'YTick',[],'NextPlot','replacechildren','Color',q(1,:));
a(3) = axes('Tag','axes3','Position',[.51 .7 .23 .25],'XTick',[],'YTick',[],'NextPlot','replacechildren','Color',q(1,:));
a(4) = axes('Tag','axes4','Position',[.76 .7 .23 .25],'XTick',[],'YTick',[],'NextPlot','replacechildren','Color',q(1,:));
a(5) = axes('Tag','axes5','Units','normalized','Position',[.01 .01 .5 .5],'XTick',[],'YTick',[],'Color',[0 0 0]);
set(a(5),'Units','pixels');   

b(1) = uicontrol('Tag','okbutton','Style','pushbutton','Units','normalized','Position',[.88 .03 .10 .04],'String','Ok','Callback','eventviewer(''okbutton_Callback'',gcbo,guidata(gcbo))');
b(2) = uicontrol('Tag','nextbutton','Style','pushbutton','Units','normalized','Position',[.85 .6 .10 .04],'String','Next','Callback','eventviewer(''nextbutton_Callback'',gcbo,guidata(gcbo))');
b(3) = uicontrol('Tag','backbutton','Style','pushbutton','Units','normalized','Position',[.75 .6 .10 .04],'String','Back','Callback','eventviewer(''backbutton_Callback'',gcbo,guidata(gcbo))');
b(4) = uicontrol('Tag','excludebutton','Style','pushbutton','Units','normalized','Position',[.60 .6 .10 .04],'String','Exclude','Callback','eventviewer(''excludebutton_Callback'',gcbo,guidata(gcbo))');
b(5) = uicontrol('Tag','unexcludebutton','Style','pushbutton','Units','normalized','Position',[.45 .6 .15 .04],'String','UnExclude','Callback','eventviewer(''unexcludebutton_Callback'',gcbo,guidata(gcbo))');
b(6) = uicontrol('Tag','toendbutton','Style','pushbutton','Units','normalized','Position',[.85 .55 .10 .04],'String','Last','Callback','eventviewer(''endbutton_Callback'',gcbo,guidata(gcbo))');
b(7) = uicontrol('Tag','tobeginningbutton','Style','pushbutton','Units','normalized','Position',[.75 .55 .10 .04],'String','First','Callback','eventviewer(''beginningbutton_Callback'',gcbo,guidata(gcbo))');
b(8) = uicontrol('Tag','viewallbutton','Style','pushbutton','Units','normalized','Position',[.33 .57 .09 .04],'String','No box','Callback','eventviewer(''viewallbutton_Callback'',gcbo,guidata(gcbo))');
b(9) = uicontrol('Tag','decreaseXbutton','Style','pushbutton','Units','normalized','Position',[.065 .615 .05 .02],'String','<','Callback','eventviewer(''decreaseXbutton_Callback'',gcbo,guidata(gcbo))');
b(10) = uicontrol('Tag','increaseXbutton','Style','pushbutton','Units','normalized','Position',[.115 .615 .05 .02],'String','>','Callback','eventviewer(''increaseXbutton_Callback'',gcbo,guidata(gcbo))');
b(11) = uicontrol('Tag','decreaseYbutton','Style','pushbutton','Units','normalized','Position',[.22 .615 .05 .02],'String','<','Callback','eventviewer(''decreaseYbutton_Callback'',gcbo,guidata(gcbo))');
b(12) = uicontrol('Tag','increaseYbutton','Style','pushbutton','Units','normalized','Position',[.27 .615 .05 .02],'String','>','Callback','eventviewer(''increaseYbutton_Callback'',gcbo,guidata(gcbo))');
b(13) = uicontrol('Tag','cancelbutton','Style','pushbutton','Units','normalized','Position',[.775 .03 .10 .04],'String','Cancel','Callback','eventviewer(''cancelbutton_Callback'',gcbo,guidata(gcbo))');

%exclude selection button not yet implemented
%b(14) = uicontrol('Tag','excludeselectedbutton','Style','pushbutton','Units','normalized','Position',[.82 .45 .15 .04],'String','Excl selection','Callback','eventviewer(''excludeselectedbutton_Callback'',gcbo,guidata(gcbo))');

e(1) = uicontrol('Tag','indexedit','Style','edit','Units','normalized','BackgroundColor',[1 1 1],'Position',[.6 .55 .10 .04],'Callback','eventviewer(''indexedit_Callback'',gcbo,guidata(gcbo))');
e(2) = uicontrol('Tag','xShowEdit','Style','edit','Units','normalized','BackgroundColor',[1 1 1],'Position',[.065 .57 .10 .04],'Callback','eventviewer(''xyShowEdit_Callback'',gcbo,guidata(gcbo))');
e(3) = uicontrol('Tag','yShowEdit','Style','edit','Units','normalized','BackgroundColor',[1 1 1],'Position',[.22 .57 .10 .04],'Callback','eventviewer(''xyShowEdit_Callback'',gcbo,guidata(gcbo))');

t(1) = uicontrol('Tag','text1','Style','text','Units','normalized','Position',[.04 .585 .02 .02],'String','X');
t(2) = uicontrol('Tag','text2','Style','text','Units','normalized','Position',[.17 .585 .05 .02],'String','Y');
t(3) = uicontrol('Tag','text3','Style','text','Units','normalized','Position',[.51 .535 .07 .05],'String','Index:');


% load the 'waves' variable that contains a SPIKELENGTH by 4 by N matrix (can use int16 to save memory)
load(clustattrib.datafile);
try
    [datasize1, datasize2, datasize3] = size(waves);
catch
    error(['The file ', clustattrib.datafile, ' contains no variable named ''waves''.']);
end
if (datasize2 ~= 4)
    error(['The ''waves'' variable in the file ',clustattrib.datafile,' must have a 2nd dimension size of 4']);
end

%the indeces to the cluster's spikes should also include previously
%excluded spikes
index = clustattrib.clusters{clustnum}.index;
try
    index = unique([index;clustattrib.eventeditindex{clustnum}(:,1)]);
end

%plot the 2 dimensional paramter space
currx = get(mainfighandles.listbox1,'Value');
curry = get(mainfighandles.listbox2,'Value');
set(e(2),'String',num2str(currx));
set(e(3),'String',num2str(curry));
axes(a(5));
hold on;
%create the currently showing points
paramplot = plot(single(clustdata.params(index,currx)),single(clustdata.params(index,curry)),'Marker','.','MarkerSize',1,'LineStyle','none','Color',[.5 .5 .5]);
%create the 'not showing' and 'excluded' points (colored white)
paramplot2 = plot(single(clustdata.params(index,currx)),single(clustdata.params(index,curry)),'Marker','.','MarkerSize',1,'LineStyle','none','Color',clustercolor);
paramplot3 = plot(single(clustdata.params(index(1),currx)),single(clustdata.params(index(1),curry)),'Marker','.','MarkerSize',1,'LineStyle','none','Color',[1 1 1]);
try
    excludes = clustattrib.eventeditindex{clustnum};
    excludeindex = excludes(:,1);
catch
    excludes = [];
    excludeindex = [];
end
set(paramplot3,'XData',single(clustdata.params(excludeindex,currx)));
set(paramplot3,'YData',single(clustdata.params(excludeindex,curry)));

%find the maximum and minimum amplitudes across all spikes (to set the
%axes)
m = max(max(clustdata.params(index,2:5)));
n = -m;
m = double(max(max(max(waves(:,:,index)),[],3)));
n = double(min(min(min(waves(:,:,index)),[],3)));
%m = m+ (.1*(m-n));
%n = n+ (.1*(m-n));
currindex = 1;

%in the parameter plot, show which spike is currently displayed
spikering = plot(clustdata.params(index(currindex),currx),clustdata.params(index(currindex),curry),'Marker','o','MarkerSize',6,'LineStyle','none','Color',[1 0 0]);
set(e(1),'String',num2str(double(index(1))));
set(e(1),'Value',double(index(1)));

%create the plot objects and initialize to the first spike
for i = 1:4
    axes(a(i));
    axis([0 datasize1+1 n m]);
    
    line(1:datasize1',zeros(datasize1,1),'Color',[0 0 0],'LineWidth',1);
    hold on
    l(i) = line(1:datasize1',waves(:,i,index(1)),'Color',[0 0 0],'LineWidth',3);
    hold on
    l2(i) = line(1:datasize1',waves(:,i,index(1)),'Color',clustercolor,'LineWidth',2);
end
%if the current spike is excluded, color it red
try
	if ismember(index(1),clustattrib.eventeditindex{clustnum}(:,1));
        set(l2,'Color',[1 1 1]);
    end
end

%attach important variables to the figure handle
fighandles = guihandles(fighandle);
fighandles.mainfighandles = mainfighandles;
fighandles.thisfig = fighandle;
fighandles.data = waves(:,:,index);
fighandles.index = index;
fighandles.allpointsindex = index;
fighandles.min = n;
fighandles.max = m;
fighandles.clustnum = clustnum;
fighandles.currindex = currindex;
fighandles.axes = a;
fighandles.l = l;
fighandles.l2 = l2;
fighandles.currx = currx;
fighandles.curry = curry;
fighandles.spikering = spikering;
fighandles.paramplot = paramplot;
fighandles.paramplot2 = paramplot2;
fighandles.paramplot3 = paramplot3;
fighandles.drawsquare = 0;
fighandles.boxon = 0;
fighandles.boxlocX = 0;
fighandles.boxlocY = 0;
fighandles.clustercolor = clustercolor;
fighandles.excludes = excludes;
fighandles.toolselect = 1;
guidata(fighandle,fighandles);
%------------------------------------------------------------
function plotwaves(currindex, fighandles)
%called when a new spike should be displayed

global clustdata;
global clustattrib;

%set the spike index indicator
set(fighandles.indexedit,'String',num2str(double(fighandles.index(currindex))));
set(fighandles.indexedit,'Value',double(fighandles.index(currindex)));

%if the spike is excluded, color it red

if (~isempty(currindex))
    set(fighandles.l2,'Visible','on');
    set(fighandles.l,'Visible','on');
    set(fighandles.spikering,'Visible','on');
    try
        if ismember(fighandles.index(currindex),fighandles.excludes(:,1));
            set(fighandles.l2,'Color',[1 1 1]);
        else
            set(fighandles.l2,'Color',fighandles.clustercolor);
        end
    catch
        set(fighandles.l2,'Color',fighandles.clustercolor);
    end
    %because fighandles.data only contains the points that belong to the
    %cluster, we need to translate the real index number.
    transindex = find(fighandles.allpointsindex == fighandles.index(currindex)); 
    
    %change the x and y data for the plot objects
    for i = 1:4
        set(fighandles.l(i),'YData',fighandles.data(:,i,transindex));
        set(fighandles.l2(i),'YData',fighandles.data(:,i,transindex));
        
    end
    
    %change the ring in the parameter plot to locate the new spike
    set(fighandles.spikering,'XData',clustdata.params(fighandles.index(currindex),fighandles.currx));
    set(fighandles.spikering,'YData',clustdata.params(fighandles.index(currindex),fighandles.curry));
else
    set(fighandles.l,'Visible','off');
    set(fighandles.l2,'Visible','off');
    set(fighandles.spikering,'Visible','off');
end
%------------------------------------------------------
function plotExcludedPoints(fighandles)

global clustattrib;
global clustdata;

try
    excludes = fighandles.excludes(:,1);
catch
    excludes = [];
end
tmpindex = setdiff(fighandles.allpointsindex,setdiff(fighandles.index,excludes));
set(fighandles.paramplot2,'XData',single(clustdata.params(fighandles.index,fighandles.currx)));
set(fighandles.paramplot2,'YData',single(clustdata.params(fighandles.index,fighandles.curry)));
set(fighandles.paramplot3,'XData',single(clustdata.params(excludes,fighandles.currx)));
set(fighandles.paramplot3,'YData',single(clustdata.params(excludes,fighandles.curry)));
%------------------------------------------------------------
function endbutton_Callback(hObject,fighandles)
%go to the last spike

global clustdata;
global clustattrib;

fighandles.currindex = length(fighandles.index);
currindex = fighandles.currindex;
plotwaves(currindex,fighandles);

guidata(fighandles.thisfig,fighandles);
set(hObject,'Selected','off');
set(fighandles.thisfig,'CurrentObject',fighandles.thisfig);
figure(fighandles.thisfig);
%------------------------------------------------------------
function beginningbutton_Callback(hObject,fighandles)
%go to the first spike

global clustdata;
global clustattrib;

fighandles.currindex = 1;
currindex = fighandles.currindex;
plotwaves(currindex,fighandles);

guidata(fighandles.thisfig,fighandles);
set(hObject,'Selected','off');
set(fighandles.thisfig,'CurrentObject',fighandles.thisfig);
figure(fighandles.thisfig);
%------------------------------------------------------------
function nextbutton_Callback(hObject,fighandles)
%go to the next spike

global clustdata;
global clustattrib;

currindex = fighandles.currindex;
if (currindex < length(fighandles.index))
    currindex = currindex + 1;
end
fighandles.currindex = currindex;
plotwaves(currindex,fighandles);

guidata(fighandles.thisfig,fighandles);
set(hObject,'Selected','off');
set(fighandles.thisfig,'CurrentObject',fighandles.thisfig);
figure(fighandles.thisfig);
%--------------------------------------------------------------
function backbutton_Callback(hObject,fighandles)
%go back one spike

global clustdata;
global clustattrib;

currindex = fighandles.currindex;
if (currindex > 1)
    currindex = currindex - 1;
end
fighandles.currindex = currindex;
plotwaves(currindex,fighandles);

guidata(fighandles.thisfig,fighandles);
set(hObject,'Selected','off');
set(fighandles.thisfig,'CurrentObject',fighandles.thisfig);
figure(fighandles.thisfig);
%--------------------------------------------------------------------
function excludebutton_Callback(hObject,fighandles)
%exclude current spike

global clustdata;
global clustattrib;

try
	if ismember(fighandles.index(fighandles.currindex),fighandles.excludes(:,1));
        %this spike has been excluded before. Set the first bit of exclues(index,2) to 1. This is an ID for this program;
        fighandles.excludes(find(fighandles.excludes(:,1) == fighandles.index(fighandles.currindex)),2) = ...
            fastbitset(fighandles.excludes(find(fighandles.excludes(:,1) == fighandles.index(fighandles.currindex)),2),1,true);
	else
        %otherwise just add the spike to the exclude list with an
        %identifyer of 1
        fighandles.excludes = [fighandles.excludes;int32([fighandles.index(fighandles.currindex) 1])];
	end
catch
    fighandles.excludes = int32([fighandles.index(fighandles.currindex) 1]);
end

plotwaves(fighandles.currindex,fighandles);
plotExcludedPoints(fighandles);
guidata(fighandles.thisfig,fighandles);
set(fighandles.thisfig,'CurrentObject',fighandles.thisfig);
figure(fighandles.thisfig);
%--------------------------------------------------------------------
function unexcludebutton_Callback(hObject,fighandles)
%un-exclude the current spike

global clustdata;
global clustattrib;

try
	if ismember(fighandles.index(fighandles.currindex),fighandles.excludes(:,1));
        %remove the exclude
        fighandles.excludes = fighandles.excludes(find(fighandles.excludes(:,1) ~= fighandles.index(fighandles.currindex)),:);
	end
end

plotwaves(fighandles.currindex,fighandles);
plotExcludedPoints(fighandles);
guidata(fighandles.thisfig,fighandles);
set(fighandles.thisfig,'CurrentObject',fighandles.thisfig);
figure(fighandles.thisfig);
%-------------------------------------------
function excludeselectedbutton_Callback(hObject,fighandles)

%-------------------------------------------------------------------
function okbutton_Callback(hObject,fighandles)
%update changes in matclust and quit when 'ok' is pressed

global clustdata;
global clustattrib;
clustnum = fighandles.clustnum;
in = fastorbit(clustattrib.pointexclude,clustattrib.clusters{clustnum}.polyindex(:,4)); %which points are excluded at any polygon
in2 = fastorbit(clustattrib.pointinclude,clustattrib.clusters{clustnum}.polyindex(:,4)); %which points have been included at any polygon
in3 = false(length(clustdata.params(:,1)),1); %which points are excluded individually

try
    if ~isempty(fighandles.excludes)
        in3(fighandles.excludes(:,1)) = true; %if any points were excluded individually, it will be stored in clustattrib.eventeditindex
        clustattrib.eventeditindex{fighandles.clustnum} = fighandles.excludes;    
    end
end

clustattrib.clusters{clustnum}.index = uint32(find(~in & in2 & ~in3));
matclust('addnewstate','edit single events',fighandles.mainfighandles); 
matclust('plotgraph',fighandles.mainfighandles); 
delete(fighandles.eventviewer);
%-----------------------------------------------------------------------
function cancelbutton_Callback(hObject,fighandles)

delete(fighandles.thisfig);
%---------------------------------------------------------------------
function indexedit_Callback(hObject,fighandles)
%this is called when the index edit box is changed manually

global clustdata;
global clustattrib;

currindex = fighandles.currindex;
nextindex = get(hObject,'String');
nextindex = str2num(nextindex);
if ((~isempty(nextindex)) & (nextindex >=1))
    [minval, currindex] = min(abs(double(fighandles.index) - nextindex)); %find the closest value to the chosen index    
end

fighandles.currindex = currindex;
plotwaves(currindex,fighandles);
guidata(fighandles.thisfig,fighandles);
set(hObject,'Selected','off');
set(fighandles.thisfig,'CurrentObject',fighandles.thisfig);
figure(fighandles.thisfig);
%---------------------------------------------------
function xyShowEdit_Callback(hObject,fighandles)

global clustdata;
global clustattrib;

currx = str2num(get(fighandles.xShowEdit,'String'));
curry = str2num(get(fighandles.yShowEdit,'String'));
numparams = size(clustdata.params,2);

if ( (~isempty(currx)) & (~isempty(curry)) & (curry >= 1) & (currx >= 1) & (currx <= numparams) & (curry <= numparams) )
    if (fighandles.boxon)
        if ( (currx == fighandles.boxlocX) & (curry == fighandles.boxlocY) )
            set(fighandles.lines,'Visible','on');
        else
            set(fighandles.lines,'Visible','off');
        end        
    end
    fighandles.currx = currx;
    fighandles.curry = curry;
    axes(fighandles.axes5);
        
    set(fighandles.paramplot,'XData',single(clustdata.params(fighandles.allpointsindex,currx)));
    set(fighandles.paramplot,'YData',single(clustdata.params(fighandles.allpointsindex,curry)));
    
    plotExcludedPoints(fighandles);
    plotwaves(fighandles.currindex,fighandles);
    guidata(fighandles.thisfig,fighandles);
else
    set(fighandles.xShowEdit,'String',num2str(fighandles.currx));
    set(fighandles.yShowEdit,'String',num2str(fighandles.curry));
end

set(fighandles.thisfig,'CurrentObject',fighandles.thisfig);
figure(fighandles.thisfig);
%---------------------------------------------------------
function viewallbutton_Callback(hObject,fighandles)
% this function erases the view box

try
    relativeindex = find(fighandles.allpointsindex == fighandles.index(fighandles.currindex));
catch
    relativeindex = 1;
end
    
fighandles.index = fighandles.allpointsindex;
fighandles.boxon = 0;
try
    delete(fighandles.lines);
end
fighandles.currindex = relativeindex;
plotExcludedPoints(fighandles);
plotwaves(fighandles.currindex,fighandles);
guidata(fighandles.thisfig,fighandles);
%----------------------------------------------------------
function decreaseXbutton_Callback(hObject,fighandles)

set(fighandles.xShowEdit,'String',num2str(fighandles.currx - 1));
xyShowEdit_Callback(hObject,fighandles);
%----------------------------------------------------------
function increaseXbutton_Callback(hObject,fighandles)

set(fighandles.xShowEdit,'String',num2str(fighandles.currx + 1));
xyShowEdit_Callback(hObject,fighandles);
%----------------------------------------------------------
function decreaseYbutton_Callback(hObject,fighandles)

set(fighandles.yShowEdit,'String',num2str(fighandles.curry - 1));
xyShowEdit_Callback(hObject,fighandles);
%----------------------------------------------------------
function increaseYbutton_Callback(hObject,fighandles)

set(fighandles.yShowEdit,'String',num2str(fighandles.curry+ 1));
xyShowEdit_Callback(hObject,fighandles);
%---------------------------------------------------
function eventviewer_WindowButtonDownFcn(hObject, fighandles)
%called whenever the left mouse button in pressed within the figure 

global clustattrib;
global clustdata;

figloc = get(fighandles.axes5,'Position');
M = get(hObject,'CurrentPoint');
%if the mouse click was within the graph window

if ((M(1)>figloc(1))&&(M(1)<(figloc(1)+figloc(3)))&&(M(2)>figloc(2))&&(M(2)<(figloc(2)+figloc(4))))
    clicktype = get(hObject,'SelectionType');
    if ((strcmp(clicktype,'extend')) && (fighandles.drawsquare == 0))
        fighandles.toolselect = 2;
    end
    if (fighandles.toolselect == 1)
        if (fighandles.drawsquare == 0) %this is the first point
            
            try
                delete(fighandles.lines);
            end
            fighandles.drawsquare = 1;
            fighandles.boxlocX = fighandles.currx;
            fighandles.boxlocY = fighandles.curry;
            fighandles.boxon = 1;
            temp = get(fighandles.axes5,'CurrentPoint');    
            fighandles.vertices = repmat(temp(1,1:2),3,1);
            fighandles.lines = line(fighandles.vertices(:,1),fighandles.vertices(:,2));
            set(fighandles.lines,'Color',[1 1 1]);       
                        
        else %we are in the middle of drawing a polygon
            
            if  strcmp(clicktype,'open')   %double click ends the drawing
                fighandles = releasepolydraw(fighandles);
                
            else
                temp = get(fighandles.axes5,'CurrentPoint');
                fighandles.vertices(end-1,1:2) = temp(1,1:2);
                fighandles.vertices(end,1:2) = temp(1,1:2);
                fighandles.vertices(end+1,1:2) = fighandles.vertices(1,1:2);
                set(fighandles.lines,'XData',fighandles.vertices(:,1));
                set(fighandles.lines,'YData',fighandles.vertices(:,2));                              
                if strcmp(clicktype,'extend')
                    fighandles = releasepolydraw(fighandles);
                end
            end                               
        end   
        
    elseif (fighandles.toolselect == 2)
        try
            delete(fighandles.lines);
        end
        %the square is set up here, and made permanent once the mouse button is released  
        fighandles.drawsquare = 1;
        fighandles.boxlocX = fighandles.currx;
        fighandles.boxlocY = fighandles.curry;
        fighandles.boxon = 1;
        temp = get(fighandles.axes5,'CurrentPoint');    
        fighandles.vertices = repmat(temp(1,1:2),5,1);
        fighandles.lines = line(fighandles.vertices(1:4,1),fighandles.vertices(1:4,2));
        set(fighandles.lines,'Color',[1 1 1]);        
    end
end
guidata(fighandles.thisfig,fighandles);
%-------------------------------------------------------------
function eventviewer_WindowButtonMotionFcn(hObject, fighandles)
%called whenever the mouse moves within the figure window

global clustattrib;
global clustdata;

figloc = get(fighandles.axes5,'Position');
M = get(hObject,'CurrentPoint');
%if the mouse is within the graph window
if ((M(1)>figloc(1))&&(M(1)<(figloc(1)+figloc(3)))&&(M(2)>figloc(2))&&(M(2)<(figloc(2)+figloc(4))))    
    
    figpoint = get(fighandles.axes5,'CurrentPoint');
    figpoint = figpoint(1,1:2);
    
    s = get(fighandles.axes5,'Position');
    
    if (fighandles.toolselect == 1)
        set(hObject,'Pointer','crosshair');         
        if (fighandles.drawsquare) %polygon currently being drawn
            fighandles.vertices(end-1,:) = figpoint;
            
            tmp = fighandles.vertices;
            set(fighandles.lines,'XData',tmp(:,1));
            set(fighandles.lines,'YData',tmp(:,2));
        end
    elseif (fighandles.toolselect == 2)   
        set(hObject,'Pointer','crosshair');         
        if (fighandles.drawsquare) %square currently being drawn
            fighandles.vertices(2,1) = figpoint(1);
            fighandles.vertices(3,:) = figpoint;
            fighandles.vertices(4,2) = figpoint(2);
            set(fighandles.lines,'XData',fighandles.vertices(:,1));
            set(fighandles.lines,'YData',fighandles.vertices(:,2));
            guidata(fighandles.thisfig,fighandles);
        end     
    end
else
    set(hObject,'Pointer','arrow');
end
%-----------------------------------------------------
function eventviewer_WindowButtonUpFcn(hObject, fighandles)
%called whenever the left mouse button is released inside the figure window

global clustattrib;
global clustdata;


try
    relativeindex = find(fighandles.allpointsindex == fighandles.index(fighandles.currindex));
catch
    relativeindex = 1;
end
if ((fighandles.toolselect == 2) && (fighandles.drawsquare)) %square was being drawn
    fighandles.drawsquare = 0;
    fighandles.toolselect = 1;
    fighandles.index = fighandles.allpointsindex(find(fastinpoly(clustdata.params(fighandles.allpointsindex,fighandles.currx), ... 
                    clustdata.params(fighandles.allpointsindex,fighandles.curry),fighandles.vertices(:,1),fighandles.vertices(:,2))));
    
    plotExcludedPoints(fighandles);
    
    if (~isempty(fighandles.index))
        %if the currently viewed point is in the box, then down change the
        %current point.  Otherwise change the index to 1;
        if ismember(fighandles.allpointsindex(relativeindex),fighandles.index)
            fighandles.currindex = find(fighandles.index == fighandles.allpointsindex(relativeindex));
        else
            fighandles.currindex = 1;
        end
    else
        fighandles.currindex = [];
    end
    plotwaves(fighandles.currindex,fighandles);
    
end
guidata(fighandles.thisfig,fighandles);
%----------------------------------------------------------------------
function eventviewer_KeyPressFcn(hObject,fighandles)
%a keyboard key was pressed

global clustdata;
global clustattrib;

try
    keypress = get(hObject,'CurrentCharacter');
    numkeypress = double(keypress); %turn the character into an ASCII number
    switch numkeypress
        case 110 %n key = next spike
            nextbutton_Callback(hObject,fighandles);
        case 31 %down arrow = next spike
            nextbutton_Callback(hObject,fighandles);
        case 98 %b key = back one spike
            backbutton_Callback(hObject,fighandles);
        case 30 %up arrow = back one spike
            backbutton_Callback(hObject,fighandles);
        case 32 %space bar toggles 'exclude' for the current spike
            
            %they may be no exluded spikes for this cluster, in which case
            % excludes will be empty.
            try
                ismember(fighandles.index(fighandles.currindex),fighandles.excludes(:,1));
                indexexists = 1;
            catch
                indexexists = 0;
            end        
            if ((indexexists) & (ismember(fighandles.index(fighandles.currindex),fighandles.excludes(:,1))));
                unexcludebutton_Callback(hObject,fighandles);
            else
                excludebutton_Callback(hObject,fighandles);
            end
    end   
end
%------------------------------------------------------------
function eventviewer_CloseRequestFnc(hObject,fighandles)
%user closed window

delete(hObject);
%-----------------------------------------------------------
function fighandles = releasepolydraw(fighandles)

global clustdata;
global clustattrib;

if ((fighandles.toolselect == 1) && (fighandles.drawsquare))
    fighandles.drawsquare = 0;
    
    fighandles.vertices(end-1,1:2) = fighandles.vertices(end,1:2);
    fighandles.vertices = fighandles.vertices(1:end-1,1:2);              
    try
        relativeindex = find(fighandles.allpointsindex == fighandles.index(fighandles.currindex));
    catch
        relativeindex = 1;
    end
    fighandles.index = fighandles.allpointsindex(find(fastinpoly(clustdata.params(fighandles.allpointsindex,fighandles.currx), ... 
        clustdata.params(fighandles.allpointsindex,fighandles.curry),fighandles.vertices(:,1),fighandles.vertices(:,2))));
    
    plotExcludedPoints(fighandles);
    
    if (~isempty(fighandles.index))
        %if the currently viewed point is in the box, then down change the
        %current point.  Otherwise change the index to 1;
        if ismember(fighandles.allpointsindex(relativeindex),fighandles.index)
            fighandles.currindex = find(fighandles.index == fighandles.allpointsindex(relativeindex));
        else
            fighandles.currindex = 1;
        end
    else
        fighandles.currindex = [];
    end
    guidata(fighandles.thisfig,fighandles);
    plotwaves(fighandles.currindex,fighandles);
    
end