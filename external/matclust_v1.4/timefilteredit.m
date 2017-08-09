function timefilteredit(varargin)

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
        'CloseRequestFcn','timefilteredit(''addtimeFigure_CloseRequestFnc'',gcbo,guidata(gcbo))', ...
        'WindowStyle','modal');

currval = get(mainfighandles.timefilterList,'Value');
names = get(mainfighandles.timefilterList,'String');

filterstring = names{currval};
if (length(filterstring)>2)
    dash = strfind(filterstring,'-');
    spaces = strfind(filterstring,' ');
    lastspace = spaces(end);
    title = filterstring(4:lastspace-1);
    startnum = filterstring(lastspace+1:dash-1);
    endnum = filterstring(dash+1:end);
else
    title = '';
    startnum = '';
    endnum = '';
end

e(1) = uicontrol('Tag','edit1','Style','edit','Position',[70 152 115 20],'BackgroundColor',[1 1 1],'String',title);
e(2) = uicontrol('Tag','edit2','Style','edit','Position',[70 102 115 20],'BackgroundColor',[1 1 1],'String',startnum);
e(3) = uicontrol('Tag','edit3','Style','edit','Position',[70 52 115 20],'BackgroundColor',[1 1 1],'String',endnum);
b(1) = uicontrol('Tag','cancelbutton','Style','pushbutton','Position',[96 5 50 20],'String','Cancel','Callback','timefilteredit(''cancelbutton_Callback'',gcbo,guidata(gcbo))');
b(2) = uicontrol('Tag','okbutton','Style','pushbutton','Position',[147 5 50 20],'String','Ok','Callback','timefilteredit(''okbutton_Callback'',gcbo,guidata(gcbo))');
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

currval = get(fighandles.mainfighandles.timefilterList,'Value');
names = get(fighandles.mainfighandles.timefilterList,'String');

filterstring = names{currval};
filterstring = filterstring(1:3);
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
    newname = [filterstring,newname];
    
    names{currval} = newname;
    %clustdata.timefiltermemmap(find(clustdata.timefiltermemmap>currval)) = clustdata.timefiltermemmap(find(clustdata.timefiltermemmap>currval))+1;
        
    clustdata.timefilterranges(currval,1:2) = [t1 t2];
    set(fighandles.mainfighandles.timefilterList,'String',names);
    clustdata.timefilternames = names;
    passfilter = ((clustdata.params(:,1) >= t1)&(clustdata.params(:,1) <= t2));
    memval = find(clustdata.timefiltermemmap==currval);
    
    clustdata.timefilters = fastbitset(clustdata.timefilters,memval,passfilter);
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

recalcpoints(currval);

matclust('addnewstate','edit time filter',fighandles.mainfighandles);

tmpfilter = fastandbit(clustdata.timefilters,find(clustdata.timefiltersOn));
tmpfilter2 = fastandbit(clustdata.otherfilters,find(clustdata.otherfiltersOn));

clustdata.filteredpoints = tmpfilter & tmpfilter2;
matclust('plotgraph',fighandles.mainfighandles);
    

delete(fighandles.addtimeFigure);
%-------------------------------------------------

function recalcpoints(source)
% allows the user to copy all polygons for cluster clustnum from time
% filter source to time filter dest.

global clustattrib;
global graphattrib;
global figattrib;
global clustdata;

%t = clustdata.timefiltersOn(:)';

set(figattrib.handles.statusText,'String','Finding points inside shape...');
set(figattrib.handles.statusText,'UserData','Finding points inside shape...');




%in = clustattrib.cluster0;  %all points allowed in current filter(s)

for clustnum = clustattrib.clustersOn'

    def = clustattrib.clusters{clustnum}.defineaxes; 


    if ~isempty(def)
    %clustattrib.clusters{dest}.defineaxes = def;
    %set(figattrib.clustcontrol(dest,3:end),'Visible','on'); %turn on cluster control    
	%if ~sum(ismember(clustattrib.clustersOn,dest))    
    %     clustattrib.clustersOn = [clustattrib.clustersOn;dest]; %add to 'on' list
    %end
        for i = 1:size(def,1)
        
            a1 = def(i,1);
            a2 = def(i,2);
            findex = def(i,3);
            sourcefilt = clustattrib.filterindex{findex};
            %destfilt = [dest sourcefilt(2:end)];
            if (sourcefilt(1) == source)

               
    
                filterval = find(clustdata.timefiltermemmap == source);
                tmpfilter = fastandbit(clustdata.timefilters,filterval);
                tmpfilter2 = fastandbit(clustdata.otherfilters,sourcefilt(2:end));
                
                in = find(tmpfilter & tmpfilter2);
                polynum = clustattrib.clusters{clustnum}.polyindex( ... 
                   find((clustattrib.clusters{clustnum}.polyindex(:,1)==a1)&(clustattrib.clusters{clustnum}.polyindex(:,2)==a2)&(clustattrib.clusters{clustnum}.polyindex(:,3)==findex)),4);

   
                memvector = ceil(polynum/32);  
                

                xpoints = ((clustdata.params(:,a1)-clustdata.datarange(1,a1))+1)/(clustdata.datarange(2,a1)-clustdata.datarange(1,a1));
                ypoints = ((clustdata.params(:,a2)-clustdata.datarange(1,a2))+1)/(clustdata.datarange(2,a2)-clustdata.datarange(1,a2));
                vertices = graphattrib.polyg(a1,a2,findex).relativevertices{clustnum};
                xvert = vertices(:,1);
                yvert = vertices(:,2);
                excluding = false(length(clustdata.params(:,1)),1);  
                including = false(length(clustdata.params(:,1)),1); %default assumption: all points are outside the polygon
                excluding(in) = true;  %default assumption: all points in the current filter are excluded by this polygon (note: points not in the filter are not excluded)
                %in = in(find(inpoly(xpoints(in),ypoints(in),xvert,yvert))); %find the points of the current filter that are inside the polygon
                in = in(find(fastinpoly(xpoints(in),ypoints(in),xvert,yvert))); %find the points of the current filter that are inside the polygon

                excluding(in) = false;  %points within polygon are not excluded
                including(in) = true; %points within polygon are included
                if (size(clustattrib.pointexclude,2)<memvector) %we have to add a new column to pointexclude and pointinclude for every 32 polygons 
                    clustattrib.pointexclude(:,end+1) = 0;
                    clustattrib.pointinclude(:,end+1) = 0;
                end
                paramlength = length(clustdata.params(:,1));
                %clustattrib.pointexclude(1:paramlength,memvector) = bitset(clustattrib.pointexclude(1:paramlength,memvector),polynum-(32*(memvector-1)),excluding);  %save the bit-wise info
                %clustattrib.pointinclude(1:paramlength,memvector) = bitset(clustattrib.pointinclude(1:paramlength,memvector),polynum-(32*(memvector-1)),including);


                clustattrib.pointexclude(1:paramlength,memvector) = fastbitset(clustattrib.pointexclude(1:paramlength,memvector),polynum-(32*(memvector-1)),excluding);  %save the bit-wise info
                clustattrib.pointinclude(1:paramlength,memvector) = fastbitset(clustattrib.pointinclude(1:paramlength,memvector),polynum-(32*(memvector-1)),including);

               
                in = fastorbit(clustattrib.pointexclude,clustattrib.clusters{clustnum}.polyindex(:,4)); %which points are excluded at any polygon
                in2 = fastorbit(clustattrib.pointinclude,clustattrib.clusters{clustnum}.polyindex(:,4)); %which points have been included at any polygon
                in3 = false(length(clustdata.params(:,1)),1); %which points are excluded individually

                try
              
                    if ~isempty(clustattrib.eventeditindex{clustnum})
                        in3(clustattrib.eventeditindex{clustnum}(:,1)) = true; %if any points were excluded individually, it will be stored in clustattrib.eventeditindex
                    end
                end
               
                clustattrib.clusters{clustnum}.index = uint32(find(~in & in2 & ~in3));
                %clustattrib.clusters{clustnum}.index = find(~in);
                %set(figattrib.clustcontrol(clustnum,2),'TooltipString',['Cluster ',num2str(clustnum),': ',num2str(length(clustattrib.clusters{clustnum}.index)),' points']);

            end
        end
    end
end
   


