function varargout = rotationfunc(varargin)


if ~ischar(varargin{1})  %launch
    fighandle = figure('MenuBar','none','Color',get(0,'DefaultUicontrolBackgroundColor'),'Position',[400 400 600 600], ... 
        'Tag','rotationviewer','NumberTitle','off','Resize','on','Name','Choose rotation...' , ... 
        'CloseRequestFcn','rotationfunc(''closeRequestFcn'',gcbo,guidata(gcbo))', ...
        'WindowStyle','normal');
    figure_fill(fighandle,varargin{1});
    %output the figure handle if wanted
	if nargout > 0
		varargout{1} = fighandle;
	end
elseif ischar(varargin{1}) % call the desired function
    
        
     [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
        
                       
end

%----------------------------------------------------------------

function figure_fill(fighandle, matclusthandles)

global clustdata;
global clustattrib;
global graphattrib;


  

set(fighandle,'WindowButtonMotionFcn','rotationfunc(''WindowButtonMotionFcn'',gcbo,guidata(gcbo))');


clustExcludePoints = false(size(clustdata.timefilters,1),1);
for i = 1:length(graphattrib.blackedOutClusters)
    clustExcludePoints(clustattrib.clusters{graphattrib.blackedOutClusters(i)}.index) = true;
end

%Any clusters currently filtered out are not displayed, also the
%current time and other filters are used

index = find(clustdata.filteredpoints & ~clustExcludePoints);
currentaxes = [2 3 4];

a(1) = axes('Tag','axes1','Units','normalized','Position',[.1 .2 .8 .8],'XTick',[],'YTick',[],'NextPlot','replacechildren','Color',[0 0 0]);
t(1) = uicontrol('Tag','text1','Style','text','Units','pixels','FontUnits','pixels','FontSize',12,'Position',[100 60 100 20],'String','Plot axes','HorizontalAlignment','center');

axes(a(1));

target = [0 0 0];
T = viewmtx(0,0,0,target);
tmpdata = clustdata.params(:,currentaxes);
tmpdata(:,4) = 1;
projection = (T*(tmpdata'))';
plothandle = plot(projection(index,1),projection(index,2),'Marker','.','MarkerSize',2,'LineStyle','none','Color',[1 1 1]);
projection = projection(:,1:2);

pan on;
zoom(1);

%create controls
b(1) = uicontrol('Tag','okbutton','Style','pushbutton','Units','pixels','Position',[40 45 50 30],'String','Add','Callback','rotationfunc(''okbutton_Callback'',gcbo, guidata(gcbo))');
b(2) = uicontrol('Tag','cancelbutton','Style','pushbutton','Units','pixels','Position',[40 10 50 30],'String','Close','Callback','rotationfunc(''closeRequestFcn'',gcbo, guidata(gcbo))');
b(3) = uicontrol('Tag','filterbutton','Style','pushbutton','Units','pixels','Position',[100 10 100 30],'String','Filter points','Callback','rotationfunc(''filter_Callback'',gcbo, guidata(gcbo))');
b(4) = uicontrol('Tag','zoominbutton','Style','pushbutton','Units','pixels','Position',[480 45 50 30],'String','+','Callback','rotationfunc(''zoomin_Callback'',gcbo, guidata(gcbo))');
b(5) = uicontrol('Tag','zoomoutbutton','Style','pushbutton','Units','pixels','Position',[480 10 50 30],'String','-','Callback','rotationfunc(''zoomout_Callback'',gcbo, guidata(gcbo))');

e(1) = uicontrol('Tag','axesedit','Style','edit','Units','pixels','Position',[100 45 100 20],'String',num2str(currentaxes),'Callback','rotationfunc(''axesedit_Callback'',gcbo, guidata(gcbo))');

s(1) = uicontrol('Tag','azslider','Style','slider','Units','pixels','Position',[220 30 200 20],'Min',-180,'Max',180,'Callback','rotationfunc(''rotateslider_Callback'',gcbo, guidata(gcbo))');
s(2) = uicontrol('Tag','elslider','Style','slider','Units','pixels','Position',[440 5 20 80],'Min',-90,'Max',90,'Callback','rotationfunc(''rotateslider_Callback'',gcbo, guidata(gcbo))');


handles = guihandles(fighandle);
handles.mainfig = fighandle;
handles.currentaxes = currentaxes;
handles.plothandle = plothandle;
handles.matclusthandles = matclusthandles;
handles.target = target;
handles.index = index;

guidata(fighandle, handles);
%-------------------------------------------------------------------------------    
function zoomin_Callback(hObject,handles)

zoom(1.25);
%---------------------------------------------------------------------------------
function zoomout_Callback(hObject,handles)

zoom(.8);
%-----------------------------------------------------------------------------

function WindowButtonMotionFcn(hObject,handles)



global clustdata;
global clustattrib;
global graphattrib;

M = get(hObject,'CurrentPoint');

if (M(2) > 80)
    
    az = get(handles.azslider,'Value');
    el = get(handles.elslider,'Value');
    target = handles.target;
    T = viewmtx(az,el,0,target);
    
    
    tmpdata = clustdata.params(handles.index,handles.currentaxes);
    tmpdata(:,4) = 1;
    
    projection = (T*(tmpdata'))';
    projection = projection(:,1:2);
    
    a = axis;
    xmidpoint = a(1) + ((a(2)-a(1))/2);
    ymidpoint = a(3) + ((a(4)-a(3))/2);
    
    D = distcalc([xmidpoint ymidpoint],projection);
    [middist, minind] = min(D);
   
    handles.target = clustdata.params(minind,handles.currentaxes);
    guidata(handles.mainfig, handles);
else
    
end
%------------------------------------------------------

function rotateslider_Callback(hObject,handles)

global clustdata;
global clustattrib;
global graphattrib;



az = get(handles.azslider,'Value');
el = get(handles.elslider,'Value');
target = handles.target;
T = viewmtx(az,el,0,target);


tmpdata = clustdata.params(handles.index,handles.currentaxes);
tmpdata(:,4) = 1;

projection = (T*(tmpdata'))';




projection = projection(:,1:2);
set(handles.plothandle,'XData',projection(:,1));
set(handles.plothandle,'YData',projection(:,2));


%-----------------------------------------------------------------------


function filter_Callback(hObject,handles)

global graphattrib;
global clustattrib;
global clustdata;

clustExcludePoints = false(size(clustdata.timefilters,1),1);
for i = 1:length(graphattrib.blackedOutClusters)
    clustExcludePoints(clustattrib.clusters{graphattrib.blackedOutClusters(i)}.index) = true;
end

%Any clusters currently filtered out are not displayed, also the
%current time and other filters are used
handles.index = find(clustdata.filteredpoints & ~clustExcludePoints);

az = get(handles.azslider,'Value');
el = get(handles.elslider,'Value');
target = handles.target;
T = viewmtx(az,el,0,target);
tmpdata = clustdata.params(handles.index,handles.currentaxes);
tmpdata(:,4) = 1;

projection = (T*(tmpdata'))';
projection = projection(:,1:2);

set(handles.plothandle,'XData',projection(:,1));
set(handles.plothandle,'YData',projection(:,2));

guidata(handles.mainfig,handles);   

%-------------------------------------------------------

function axesedit_Callback(hObject,handles)

global clustdata;
global graphattrib;
global clustattrib;

clustExcludePoints = false(size(clustdata.timefilters,1),1);
for i = 1:length(graphattrib.blackedOutClusters)
    clustExcludePoints(clustattrib.clusters{graphattrib.blackedOutClusters(i)}.index) = true;
end

%Any clusters currently filtered out are not displayed, also the
%current time and other filters are used
index = find(clustdata.filteredpoints & ~clustExcludePoints);

oldaxes = handles.currentaxes;
maxparam = size(clustdata.params,2);
try
    newaxes = str2num(get(hObject,'String'));
    if (length(newaxes) ~= 3)
        error
    end
    test1 = sum(newaxes > maxparam);
    test2 = sum(newaxes < 1);
    if (test1 || test2)
        error
    end
    handles.currentaxes = newaxes;
    
    az = get(handles.azslider,'Value');
    el = get(handles.elslider,'Value');
    target = handles.target;
    T = viewmtx(az,el,0,target);
    tmpdata = clustdata.params(handles.index,handles.currentaxes);
    tmpdata(:,4) = 1;
    
    projection = (T*(tmpdata'))';
    projection = projection(:,1:2);
    
    set(handles.plothandle,'XData',projection(:,1));
    set(handles.plothandle,'YData',projection(:,2));
       
    guidata(handles.mainfig,handles);
catch
   
    set(hObject,'String',num2str(oldaxes));
end
%--------------------------------------------------------------
function okbutton_Callback(hObject,handles)

global clustdata;
global graphattrib;

% T = view(handles.axes1);
% [az, el] = view(handles.axes1);
% viewparam = get(handles.axes1,'View');
% target = get(handles.axes1,'CameraTarget');
% viewangle = get(handles.axes1,'CameraViewAngle');
%T = viewmtx(az,el,0,target);
%T = viewmtx(az,el);


az = get(handles.azslider,'Value');
el = get(handles.elslider,'Value');
target = handles.target;
T = viewmtx(az,el,0,target);

mhandles = handles.matclusthandles;

tmpdata = clustdata.params(:,handles.currentaxes);
tmpdata(:,4) = 1;

projection = (T*(tmpdata'))';


%projection = projection(:,1:2)./repmat(projection(:,4),[1 2]);
projection = projection(:,1:2);

clear tmpdata;
S = get(mhandles.listbox1,'String');

tmpfill = [clustdata.filledparam 0];
firstzero = min(find(tmpfill == 0));
%clustdata.params(:,firstzero) = outparam(:,i);
clustdata.names{firstzero} = ['ROTATION: ',num2str(handles.currentaxes)];
clustdata.filledparam(firstzero) = 1;
clustdata.rotation(firstzero).matrix = T;
clustdata.rotation(firstzero).params = handles.currentaxes;

clustdata.rotation(firstzero).datarange(1:2,1) = [min(projection(:,1));max(projection(:,1))]; %stores the minimum and maximum values (-/+ 10% of range) of each parameter
clustdata.rotation(firstzero).datarange(1:2,1) = [min(projection(:,1))-(.1*diff(clustdata.rotation(firstzero).datarange(1:2,1)));max(projection(:,1))+(.1*diff(clustdata.rotation(firstzero).datarange(1:2,1)))];
clustdata.rotation(firstzero).datarange(1:2,2) = [min(projection(:,2));max(projection(:,2))]; %stores the minimum and maximum values (-/+ 10% of range) of each parameter
clustdata.rotation(firstzero).datarange(1:2,2) = [min(projection(:,2))-(.1*diff(clustdata.rotation(firstzero).datarange(1:2,2)));max(projection(:,2))+(.1*diff(clustdata.rotation(firstzero).datarange(1:2,2)))];

S{firstzero} = ['ROTATION: ',num2str(handles.currentaxes)];
maxfilledparam = max(find(clustdata.filledparam));
graphattrib.viewbox(firstzero,1:4,1:size(graphattrib.viewbox,3)) = repmat([0 1 0 1],[1 1 size(graphattrib.viewbox,3)]);
graphattrib.viewbox(:,1:4,firstzero) = repmat([0 1 0 1],[size(graphattrib.viewbox,1) 1 1]);


set(mhandles.listbox1,'String',S);
set(mhandles.listbox2,'String',S);

matclust('addnewstate','add rotation',mhandles); 

%---------------------------------------------------------------
function closeRequestFcn(hObject,handles)

global figattrib;


delete(handles.mainfig);
figattrib.rotationwindowHandle = [];
%--------------------------------------------------------------

function [d] = distcalc(x1, x2)

if ((size(x1,1) == 1) & (size(x2,1) > 1))
	tmpx1 = zeros(size(x2,1), size(x2,2));
	tmpx1(:,1) = x1(1);
	tmpx1(:,2) = x1(2);
	x1 = tmpx1;	
end

d = sqrt(sum(((x1 - x2).^2)')'); 
        