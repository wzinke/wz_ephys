%    |\   /|   /\  ------ .---  |   |   |  .----  ------
% __ | \_/ |  /--\    |   |     |   |   | /___..    | _______________
%    |     | /    \   |   |____ |__ |___| ____/     |
%
% By Mattias Karlsson
%
% Matclust is a GUI program that runs in Matlab and allows the user to draw 
% polygons or boxes around multidimensional data points.  This 'clustering' of data is used 
% often for tetrode recordings of neural data.
% 
% Starting Matclust
% 
% Once Matclust is installed on the computer, simply type 'matclust' in the Matlab 
% prompt to start the program.  You can also feed Matclust an input variable by typing 
% 'matclust(yourvariable)'.  The input variable can be of a few different formats:
% 
% 
% 1)	a matrix.  If the input is simply a matrix, then matclust will consider each row a data 
% point, and each column represents a different parameter to cluster.  The first column 
% must be time (in seconds) in order for the graph labels to be correct.  Matclust will 
% name each parameter Column 1, Column2, ...
% 2)	a structure.  Here, one of the fields must be called 'params', which contains the matrix 
% described above.  You can also include a few other fields: 
% 
% paramnames  -- an N by 1 cell array, with N equal to the number columns in the data 
% matrix.  Each cell contains a string with the name of the corresponding parameter.
% 
% filename - this is a string with the name of a file containing original data. User written 
% helper programs can use this info if needed.  In particular, if the user wants to be able 
% to view neural spikes using eventviewer.m, then the file must contain a variable named 
% 'waves' of size SPIKELENGTH by 4 by N, where SPIKELENGTH is the amount of 
% points collected for each spike, and N is the amount of spikes.  This N should match 
% the N in the parameter matrix.  The second dimension must be 4 because eventviewer 
% assumes you are using tetrodes.  The matrix does not have to be a double array-it can 
% be int16, for example, to save memory. Eventviewever is called from the 'tools' menu in 
% the clusters' right-click menu.
% 
% customvar - this field can contain any data that user-written helper programs can use. 
% This data will be stored in the global variable 'clustdata' as clustdata.customvar.

function varargout = matclust(varargin)
% MatClust Application M-file
%    FIG = matclust
%    matclust('callback_name', ...) invoke the named callback.

DEBUGMODE = 1;  %change to 1 to pipe errors to matlab 

if ((nargin == 0)|((nargin == 1)&&(~ischar(varargin{1}))))  %launch matclust
    global figattrib;
	%create a blank figure
    if isempty(figattrib)
        
        figureColor = get(0,'DefaultUicontrolBackgroundColor');
        figureColor = [.91 .90 .83];
         fig = figure('MenuBar','none','Color',figureColor ,'Position',[270 150 750 500], ... 
        'Tag','figure1','NumberTitle','off','Name','MATCLUST Version 1.4','Renderer','OpenGL','Visible','off');
    else
        %Multiple instances of matclust are not allowed because of global
        %variable interference.
        error('MatClust may already be open -- type CLEAR GLOBAL if not')
    end
    %fill the figure with all gui's
    try
        figure1_fill(fig);
    catch
        %if figure1_fill fails, then delete the figure, display the
        %error message, and clear the created global variables.
        disp('There was a problem starting matclust.  Error in figure1_fill.');
        disp(lasterr);
        delete(fig);
        clear global clustdata graphattrib clustattrib figattrib;
        return
    end
    if (nargin == 1)
        openinputparam_Callback(guidata(gcf),varargin{1});
    end
    
    %output the figure handle if wanted
	if nargout > 0
		varargout{1} = fig;
	end
elseif ischar(varargin{1}) % call the desired function
    if (DEBUGMODE) %errors are piped to matlab 
        [varargout{1:nargout}] = feval(varargin{:});
    else %errors are piped to gui 
        try
            [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
        catch
            
            if (toc > 1)
                tic
                matclust_error;
            else
                rethrow(lasterror);
            end

        end
    end
elseif ischar(varargin{3}) % call the desired function
    if (DEBUGMODE) %errors are piped to matlab 
        [varargout{1:nargout}] = feval(varargin{3},varargin{1:2});
    else %errors are piped to gui 
        try
            [varargout{1:nargout}] = feval(varargin{3},varargin{1:2}); % FEVAL switchyard
        catch
            
            if (toc > 1)
                tic
                matclust_error;
            else
                rethrow(lasterror);
            end

        end
    end

    
end

% Section for object create functions
%--------------------------------------------------------------------
% --------------------------------------------------------------------
function figure1_fill(h)
% fills the figure with all buttons, axes, etc and initializes all variables.

%Because matlab does not allow passing by reference, this program makes use
%of matlab's global variable system instead.  This allows access to
%large variables by all functions without copying them with each function call.
%While this saves memory, it also makes opening multiple matclust instances
%problematic.  There are four global variables used:  

global clustdata;   %stores data
global graphattrib; %stores plotting and graphics related info
global clustattrib; %stores information related to the clusters and cluster boxes 
global figattrib;   %stores general figure-related information

tic;
fprintf('loading')
fprintf('.')
%initialize the fields within the global variables.
%*******************************************************************
MainFigHandle = h; %handle to the main figure
graphattrib.Himage = 0; %handle to the image displayed in the axes
graphattrib.polydraw = 0;  %used to show that a polygon is currently being drawn
graphattrib.drawsquare = 0; %used to show that a square is currently being drawn
graphattrib.magdrawsquare = 0; %used to show that a magnifying square is currently being drawn
graphattrib.magsquare = []; %contains the handle info for the magnifying square
graphattrib.handdown = 0;%used to show that the hand tool is pressed down
graphattrib.handinfo = []; %used to calculate how much the hand has dragged the image
graphattrib.polyg = []; %contains the vertices and handles of the drawn polygons
graphattrib.currentpolyhighlight = []; %handle to the polygon that is currently highlighted
graphattrib.squarehighlighton = 0; %shows that the mouse is currently pressing on one of the vertices of a polygon
graphattrib.pointmoved = 0; %shows that the user moved one the the vertices of a polygon, and the cluster needs to be redrawn
graphattrib.plothighlightclear = 1; %decides whether the polygon highlighting should be removed when the axes are redrawn
graphattrib.polypress = 0; %shows if the mouse is currently pressing on a polygon (not one of the vertices)
graphattrib.relativepolypress = []; %used while dragging polygons- stores the vertex locations relative to the mouse 
graphattrib.shiftclick = 0;
figattrib.shiftpress = 0;
figattrib.cashedpolygon = [];
figattrib.mainFigFocus = 0;
graphattrib.rightclick = 0;
clustattrib.clusters = []; %contains the axes info and polygon vertices of the clustattrib.clusters
clustattrib.filterindex = []; %indexes the filters used for each cluster box or polygon
clustattrib.takenpolys = []; %stores which columns in the 32-bit matrix of clustattrib.pointexclude are currently being used
clustattrib.eventeditindex = []; %stores which points have been individually excluded from each cluster
figattrib.toolselect = 1; %stores which tool is currently selected
currx = 1; %sets the initial x and y axes
curry = 2;
figattrib.clustcontrol = []; %array of handles to the cluster control buttons
figattrib.tFiltAllowMultiple = 0; %whether or not to allow the user to pick multiple time filters simultaneously
figattrib.polygonsDragTogether = 1; %when multple polygons from the same cluster are displayed, do they drag together
clustattrib.clustersOn = []; %stores which clustattrib.clusters are on, and their display order 
clustdata.filledparam = [];
clustattrib.currentfilepath = [];
clustattrib.currentfilename = []; %stores the name of the current MatClust file
clustattrib.currentparamfilename = [];
clustattrib.nodata = 1; %equals 1 if no file is open
clustattrib.datafile = []; %file name of .mat file containing entire waveforms
figattrib.openfiles = [];  %list of currently open files
figattrib.currentopenfile = []; %which file out of the list is currently displayed
figattrib.switchfileMenu = []; %handles to the menus that switch the currently active matclust file
figattrib.rotationwindowHandle = [];
clustattrib.lastaction = [];  %string describing the last user action
clustattrib.newchanges = 0; %0 or 1 depending on if the newest changes have been saved
clustattrib.states = []; %contains the state history of the file
clustattrib.currstate = 1;
figattrib.panelDrag = 0; %currently dragging side panel?
%*********************************************************************

%set data information
%****************************
clustdata.params = [1 1;2 2]; %set up fake data for the blank startup screen
clustdata.origparams = 1;
clustdata.names = {'',
                    ' '};
clustdata.timefilterranges = [];
clustdata.timefilters = int32(zeros(size(clustdata.params,1),1)); %contains the filtere points of up to 32 time filters
clustdata.timefilters = fastbitset(clustdata.timefilters,1,logical(1)); %the first filter is all times
clustdata.timefiltermemmap = zeros(32,1); %because the flter order can be changed by the user, this map translates the filter numbers to those stored in memory 
clustdata.timefiltermemmap(1) = 1;
clustdata.timefiltersOn = clustdata.timefiltermemmap'; %which time filters are currently on
clustdata.otherfilters = int32(zeros(size(clustdata.params,1),1));
clustdata.otherfilters = fastbitset(clustdata.otherfilters,1,logical(1));
clustdata.otherfiltermemmap = zeros(32,1);
clustdata.otherfiltermemmap(1) = 1;
clustdata.otherfiltersOn = clustdata.otherfiltermemmap'; %for the 'other' filters, more than one filter can be on at the same time 
clustdata.filtermemmap = [clustdata.timefiltermemmap; clustdata.otherfiltermemmap];
clustdata.filteredpoints = logical(ones(size(clustdata.params,1),1)); %stores which points are currently let through the filters
clustdata.datarange = [min(clustdata.params,[],1);max(clustdata.params,[],1)]; %stores the minimum and maximum values (-/+ 10% of range) of each parameter
%clustdata.datarange = [min(clustdata.params)-(.1*diff(clustdata.datarange));max(clustdata.params)+(.1*diff(clustdata.datarange))];
graphattrib.viewbox = repmat([0 1 0 1],[size(clustdata.params,2) 1 size(clustdata.params,2)]);  %stores the current view window for each axis pair [xmin xmax ymin ymax]
graphattrib.oldviewbox = graphattrib.viewbox;
clustattrib.pointexclude = int32(zeros(size(clustdata.params,1),1)); %contains double-precision vectors storing the outside-inside status of each data point for 32 different polygons
clustattrib.pointinclude = int32(zeros(size(clustdata.params,1),1));
%******************************

%set default figure preferences
%******************************
load matclust_defaults;
clustdata.UnitsPerSec = matclust_defaults.UnitsPerSec; %controls how many units in the time column of params equals one second
clustattrib.cluster0attrib.color = matclust_defaults.Cluster0Color; %stores the color of cluster 0
graphattrib.backgroundcolor = matclust_defaults.GraphBackgroundColor; %stores the color of the graphs background
figattrib.mixcolor = matclust_defaults.ClusterColors; %stores the colors of the clusters
figattrib.maxundos = matclust_defaults.MaxUndos+1; %maximum number of undo steps stored
graphattrib.resolutionfactor = matclust_defaults.ResFactor; %controls how big the spike points are
figattrib.datafoldername = matclust_defaults.DataFolder;
figattrib.foldername = matclust_defaults.UserFolder;

try
    figattrib.position = matclust_defaults.FigPosition;
catch
    figattrib.position = [270 150 750 500];
end
try
    figattrib.sidePanelWidth = matclust_defaults.sidePanelWidth; %the width of the side control panel
catch
    figattrib.sidePanelWidth = 285;
end
try
    figattrib.filterBoxDivider = matclust_defaults.filterBoxDivider; %the division point between the two filter list boxes (range between 0 and 1) 
catch
    figattrib.filterBoxDivider = .6;
end
%******************************
fprintf('.')
%set cluster info
%******************************
%clustattrib.cluster0 = [1:size(clustdata.params,1)]'; %index to the points in cluster 0
clustattrib.cluster0attrib.show = 1; %controls whether cluster 0 is visible
clustattrib.hiddenclusters = zeros(size(figattrib.mixcolor,1),2); %stores which clusters are hidden or excluded
figattrib.colors = [];
%*******************************

%set figure specs and callbacks
%********************************
set(h,'Renderer','zbuffer');
set(h,'Position',figattrib.position);
figsize = get(h,'Position');
set(h,'ResizeFcn','matclust(''figure1_ResizeFcn'',gcbo,guidata(gcbo))');
set(h,'CloseRequestFcn','matclust(''figure1_CloseRequestFcn'',gcbo,guidata(gcbo))');
set(h,'WindowButtonDownFcn','matclust(''figure1_WindowButtonDownFcn'',gcbo,guidata(gcbo))');
set(h,'WindowButtonUpFcn','matclust(''figure1_WindowButtonUpFcn'',gcbo,guidata(gcbo))');
set(h,'WindowButtonMotionFcn','matclust(''figure1_WindowButtonMotionFcn'',gcbo,guidata(gcbo))');
set(h,'KeyPressFcn',{@matclust,'figure1_KeyPressFcn'});
set(h,'KeyReleaseFcn',{@matclust,'figure1_KeyReleaseFcn'});
set(h,'WindowScrollWheelFcn',{@matclust,'figure1_WindowScrollWheelFcn'});
figattrib.figureColor = get(h,'Color');
%********************************
fprintf('.')
%create menu items
%**********************************

%Polygon menu
polycontext = uicontextmenu('Tag','polyContext');
pcontext1 = uimenu(polycontext,'Tag','polydeleteContextMenu','Label','Delete polygon','Callback','matclust(''deletepolymenuCallback'',gcbo,guidata(gcbo))');
pcontext2 = uimenu(polycontext,'Tag','addpolyfilterContextMenu','Label','Create filter','Callback','matclust(''addPolygonFilter_Callback'',gcbo,guidata(gcbo))');

%main menu
FigMenu(1) = uimenu('Label','File','Tag','fileMenu');
  FileMenu(1) = uimenu(FigMenu(1),'Label','Open','Tag','file_openMenu');
    OpenMenu(1) = uimenu(FileMenu(1),'Label','Raw parameter file','Tag','open_rawparamMenu','Callback','matclust(''openrawparam_Callback'',gcbo,guidata(gcbo))','UserData',-1);
    OpenMenu(2) = uimenu(FileMenu(1),'Label','MatClust file','Tag','open_clustfileMenu','Callback','matclust(''openclustfile_Callback'',gcbo,guidata(gcbo))','UserData',-1);  
  FileMenu(2) = uimenu(FigMenu(1),'Label','Close','Tag','file_closeMenu','Callback','matclust(''closeMenu_Callback'',gcbo,guidata(gcbo))','Enable','off');
  FileMenu(3) = uimenu(FigMenu(1),'Label','Save','Tag','file_saveMenu','Callback','matclust(''saveMenu_Callback'',gcbo,guidata(gcbo))','Separator','on','Enable','off','Accelerator','f');
  FileMenu(4) = uimenu(FigMenu(1),'Label','Save As...','Tag','file_saveasMenu','Callback','matclust(''saveasMenu_Callback'',gcbo,guidata(gcbo))','Enable','off');
  FileMenu(5) = uimenu(FigMenu(1),'Label','Export/Import','Tag','file_exportMenu','Separator','on','Enable','off');
    ExportMenu(1) = uimenu(FileMenu(5),'Label','Export time filters','Tag','exporttimefiltMenu','Callback','matclust(''exporttimefiltMenu_Callback'',guidata(gcbo))');
    ExportMenu(2) = uimenu(FileMenu(5),'Label','Import time filters','Tag','importtimefiltMenu','Callback','matclust(''importtimefiltMenu_Callback'',guidata(gcbo))');
    ExportMenu(3) = uimenu(FileMenu(5),'Label','Export cluster data','Tag','excportclustdataMenu','Callback','matclust(''exportdataMenu_Callback'',guidata(gcbo))');
    ExportMenu(4) = uimenu(FileMenu(5),'Label','Export current image','Tag','excportimageMenu','Callback','matclust(''exportimageMenu_Callback'',guidata(gcbo))');
  FileMenu(6) = uimenu(FigMenu(1),'Label','Currently open files','Tag','file_curropenMenu','Separator','on','Enable','off');
  FileMenu(7) = uimenu(FigMenu(1),'Label','Quit','Tag','file_quitMenu','Separator','on','Callback','matclust(''figure1_CloseRequestFcn'',gcbo,guidata(gcbo))');
FigMenu(2) = uimenu('Label','Edit','Tag','editMenu');
  EditMenu(1) = uimenu(FigMenu(2),'Label','Undo','Tag','undoMenu','Callback','matclust(''undoMenu_Callback'',gcbo,guidata(gcbo))','Enable','off','Accelerator','z');
  EditMenu(2) = uimenu(FigMenu(2),'Label','Redo','Tag','redoMenu','Callback','matclust(''redoMenu_Callback'',gcbo,guidata(gcbo))','Enable','off','Accelerator','r');
  EditMedu(3) = uimenu(FigMenu(2),'Label','Polygon','Tag','mainPolygonMenu','Enable','on');
    PolyMenu(1) = uimenu(EditMedu(3),'Tag','polydeletemenu','Label','Delete polygon  (delete key)','Callback','matclust(''deletepolymenuCallback'',gcbo,guidata(gcbo))','Enable','off');
    PolyMenu(2) = uimenu(EditMedu(3),'Tag','addPolygonFilterMainmenu','Label','Create filter','Callback','matclust(''addPolygonFilter_Callback'',gcbo,guidata(gcbo))','Enable','off');
    PolyMenu(3) = uimenu(EditMedu(3),'Tag','copyPolygonFilterMainmenu','Label','Copy  (shift+hold)','Callback','matclust(''cashPolygon'',gcbo,guidata(gcbo))','Enable','off','Accelerator','c');
    PolyMenu(4) = uimenu(EditMedu(3),'Tag','pastePolygonFilterMainmenu','Label','Paste  (release shift)','Callback','matclust(''pasteCashedPolygon'',gcbo,guidata(gcbo))','Enable','off','Accelerator','v');

  EditMenu(4) = uimenu(FigMenu(2),'Label','Cluster','Tag','mainClustMenu','Enable','off');
    i = 1;
    mainClusterMenu(1) = uimenu(EditMenu(4),'Tag','clustcopyMenu','Label','Copy cluster to...','Callback','matclust(''copyclustMenuFcn'',get(gcbo,''UserData''),guidata(gcbo))','UserData',i);
    mainClusterMenu(2) = uimenu(EditMenu(4),'Tag','epochcopyMenu','Label','Copy all polygons in a time filter','Callback','matclust(''copyepochMenuFcn'',get(gcbo,''UserData''),guidata(gcbo))','UserData',i);  
    mainClusterMenu(3) = uimenu(EditMenu(4),'Tag','clustdeleteMenu','Label','Delete','UserData',i);
    mainClusterMenu(4) = uimenu(mainClusterMenu(3),'Tag','clustdelallMenu','Label','Entire cluster','Callback','matclust(''clustdelete'',get(gcbo,''UserData''),guidata(gcbo),1)','UserData',i);
    mainClusterMenu(5) = uimenu(mainClusterMenu(3),'Tag','clustdelinfiltMenu','Label','Only in current time filter(s)','Callback','matclust(''clustfiltdelete'',get(gcbo,''UserData''),guidata(gcbo))','UserData',i);
    mainClusterMenu(6) = uimenu(mainClusterMenu(3),'Tag','clustdelwindowtMenu','Label','Polygon(s) in current window','Callback','matclust(''clustwindelete'',get(gcbo,''UserData''),guidata(gcbo))','UserData',i);
    mainClusterMenu(7) = uimenu(mainClusterMenu(3),'Tag','clustdelsingexcludesMenu','Label','Single point excludes','Callback','matclust(''spexcludedelete'',get(gcbo,''UserData''),guidata(gcbo))','UserData',i);
    mainClusterMenu(8) = uimenu(EditMenu(4),'Tag','colorMenu','Label','Change color','Callback','getcolor(gcbo,guidata(gcbo))','UserData',i);
    mainClusterMenu(9) = uimenu(EditMenu(4),'Tag','orderMenu','Label','Order','UserData',i);
    mainClusterMenu(10) = uimenu(mainClusterMenu(9),'Label','To Front','UserData',i,'Callback','matclust(''tofrontFcn'',gcbo,guidata(gcbo))');
    mainClusterMenu(11) = uimenu(mainClusterMenu(9),'Label','To Back','UserData',i,'Callback','matclust(''tobackFcn'',gcbo,guidata(gcbo))');
    mainClusterMenu(12) = uimenu(mainClusterMenu(9),'Label','Foreward One','UserData',i,'Callback','matclust(''forewardoneFcn'',gcbo,guidata(gcbo))');
    mainClusterMenu(13) = uimenu(mainClusterMenu(9),'Label','Back One','UserData',i,'Callback','matclust(''backoneFcn'',gcbo,guidata(gcbo))');
    mainClusterMenu(14) = uimenu(EditMenu(4),'Tag','ClustToolsMenu','Label','Tools','UserData',i);
    clusttoolsmenuNum = 14;
    mainClusterMenu(15) = uimenu(mainClusterMenu(clusttoolsmenuNum),'Tag','viewISIMenu','Label','View ISI','Callback','ISIviewer(get(gcbo,''UserData''),guidata(gcbo))','UserData',i);
    mainClusterMenu(16) = uimenu(mainClusterMenu(clusttoolsmenuNum),'Tag','viewEvents','Label','View spikes','Callback','eventviewer(get(gcbo,''UserData''),guidata(gcbo))','UserData',i);
      
    % add any additional tool functions from the 'ClustTools' folder
	currdir = pwd;
	cd(figattrib.foldername);
    cd('ClustTools');
	toolsdir = dir;
	numExtraTools = 0;
    for j = 3:length(toolsdir)
        f = findstr(toolsdir(j).name,'.');
        tmpname = toolsdir(j).name(1:f-1);
        numExtraTools = numExtraTools+1;
        mainClusterMenu(clusttoolsmenuNum+2+numExtraTools) = uimenu(mainClusterMenu(clusttoolsmenuNum),'Label',tmpname,'Callback',['matclust(''ExecuteClustTool'',''',tmpname,''',get(gcbo,''UserData''))'],'UserData',i);
	end
	cd(currdir);
    
    
FigMenu(3) = uimenu('Label','View','Tag','viewMenu');
  ViewMenu(1) = uimenu(FigMenu(3),'Label','Background color','Tag','view_bcolorMenu','Callback','getcolor(gcbo,guidata(gcbo))','UserData',-1);
  ViewMenu(2) = uimenu(FigMenu(3),'Label','Point Size','Tag','view_psizeMenu');
    PointSizeMenu(1) = uimenu(ViewMenu(2),'Label','Small','Tag','smallpsizeMenu','Callback','matclust(''psizeMenuFcn'',gcbo,guidata(gcbo))','UserData',1,'Checked','on','Accelerator','p');
    PointSizeMenu(2) = uimenu(ViewMenu(2),'Label','Large','Tag','largepsizeMenu','Callback','matclust(''psizeMenuFcn'',gcbo,guidata(gcbo))','UserData',2,'Accelerator','l');
  ViewMenu(3) = uimenu(FigMenu(3),'Label','Zoom','Tag','view_zoomMenu');
    ZoomMenu(1) = uimenu(ViewMenu(3),'Label','Zoom full','Tag','view_zoomoutMenu','Callback','matclust(''zoomFull'')');
    ZoomMenu(2) = uimenu(ViewMenu(3),'Label','Zoom all in','Tag','view_zoomallinMenu','Callback','matclust(''zoomAllIn'')');
    ZoomMenu(3) = uimenu(ViewMenu(3),'Label','Zoom all out','Tag','view_zoomalloutMenu','Callback','matclust(''zoomAllOut'')');
    ZoomMenu(4) = uimenu(ViewMenu(3),'Label','Zoom all to data','Tag','view_zoomalldataMenu','Callback','matclust(''setZoomRange'',99.9);matclust(''zoomAllOut'')');
  ViewMenu(4) = uimenu(FigMenu(3),'Label','Filter Boxes','Tag','view_fboxMenu');
    fboxMenu(1) = uimenu(ViewMenu(4),'Label','Time box bigger','Tag','view_tboxbigMenu','Callback','matclust(''setFilterBoxDivider'',.6)');
    fboxMenu(2) = uimenu(ViewMenu(4),'Label','Equal size','Tag','view_equalsizeMenu','Callback','matclust(''setFilterBoxDivider'',.5)');
    fboxMenu(3) = uimenu(ViewMenu(4),'Label','Time box smaller','Tag','view_tboxsmallMenu','Callback','matclust(''setFilterBoxDivider'',.4)');
FigMenu(4) = uimenu('Label','Tools','Tag','toolsMenu');
  ToolsMenu(1) = uimenu(FigMenu(4),'Label','Control mode','Tag','tools_controlModeMenu');
    cModeMenu(1) = uimenu(ToolsMenu(1),'Label','Draw polygon','Tag','cMode_polygonMenu','Callback','matclust(''polygonbutton_Callback'')','Accelerator','1');
    cModeMenu(2) = uimenu(ToolsMenu(1),'Label','Draw square','Tag','cMode_polygonMenu','Callback','matclust(''squarebutton_Callback'')','Accelerator','2');
    cModeMenu(3) = uimenu(ToolsMenu(1),'Label','Selection tool','Tag','cMode_polygonMenu','Callback','matclust(''arrowbutton_Callback'')','Accelerator','3');
    cModeMenu(4) = uimenu(ToolsMenu(1),'Label','Magnifier','Tag','cMode_polygonMenu','Callback','matclust(''magbutton_Callback'')','Accelerator','4');
    cModeMenu(5) = uimenu(ToolsMenu(1),'Label','Hand tool','Tag','cMode_polygonMenu','Callback','matclust(''handbutton_Callback'')','Accelerator','5');
    cModeMenu(6) = uimenu(ToolsMenu(1),'Label','Magic wand','Tag','cMode_wandMenu','Callback','matclust(''wandbutton_Callback'')','Accelerator','6');
  ToolsMenu(2) = uimenu(FigMenu(4),'Label','Analyze overlap','Tag','overlapMenu','Callback','analyzeoverlap(1);');
FigMenu(5) = uimenu('Label','Options','Tag','toolsMenu');
  OptionsMenu(1) = uimenu(FigMenu(5),'Label','Time selection behavior','Tag','options_tfiltMenu');
      TFiltMenu(1) = uimenu(OptionsMenu(1),'Label','One at a time only','Tag','oneTFiltMenu','Callback','matclust(''TFiltBehaveFcn'',gcbo,guidata(gcbo))','UserData',1,'Checked','on');
      TFiltMenu(2) = uimenu(OptionsMenu(1),'Label','Allow multiple activations','Tag','multTFiltMenu','Callback','matclust(''TFiltBehaveFcn'',gcbo,guidata(gcbo))','UserData',2,'Checked','off');
  OptionsMenu(2) = uimenu(FigMenu(5),'Label','Time units','Tag','options_tunitsMenu','Callback','matclust(''TUnitsMenuFcn'',gcbo,guidata(gcbo))');
FigMenu(6) = uimenu('Label','Help','Tag','helpMenu');
    HelpMenu(1) = uimenu(FigMenu(6),'Label','MATCLUST help','Tag','matclusthelpMenu','Callback','matclust(''matclusthelpMenuFcn'',guidata(gcbo))');
% add any additional tool functions from the 'Tools' folder
currdir = pwd;
cd(figattrib.foldername);
cd('Tools');
toolsdir = dir;
for i = 3:length(toolsdir)
    f = findstr(toolsdir(i).name,'.');
    tmpname = toolsdir(i).name(1:f-1);
    uimenu(FigMenu(4),'Label',tmpname,'Callback',['matclust(''ExecuteTool'',''',tmpname,''')']);
end
cd(currdir);
%**********************************
fprintf('.')
%create figattrib.colors (not yet implemented for later use)
%**********************************
figattrib.mixcolor2 = [1 1 1;figattrib.mixcolor];
for j = 1:length(figattrib.mixcolor2)
	for i = 1:65
        tempfigattrib.colors(i,1:3) = (figattrib.mixcolor2(j,:)*(65-i)+((figattrib.mixcolor2(j,:)+2*[.5 0 0])/3)*(i))/65;
	end
	tempfigattrib.colors = [0 0 0;tempfigattrib.colors];
    figattrib.colors = [figattrib.colors;tempfigattrib.colors];
end
%**********************************

%create the axes and toolbar buttons
%********************************************
%create graph axes
graphattrib.graphwindow = axes('Tag','axes3','Units','pixels','Position',[40 35 figsize(3)-310 figsize(4)-56],'Color',graphattrib.backgroundcolor);
buttonsep = 28;

%bcolor = get(0,'DefaultUicontrolBackgroundColor')*255;
bcolor = figattrib.figureColor*255;
bcolor2 = [200 200 200];
%create toolbar buttons
temppic = imread('poly','bmp');
temppic = changepicbackground(temppic,bcolor2);
toolbutton(1) = uicontrol('Style','togglebutton','Tag','polygonbutton','Position',[2+(buttonsep*0) figsize(4)-20 30 20],...
        'Callback','matclust(''polygonbutton_Callback'')','CData',temppic,'TooltipString','Draw polygon  Ctrl+1');
temppic = imread('square','bmp');
temppic = changepicbackground(temppic,bcolor);
toolbutton(2) = uicontrol('Style','togglebutton','Tag','squarebutton','Position',[2+(buttonsep*1) figsize(4)-20 30 20],...
        'Callback','matclust(''squarebutton_Callback'')','CData',temppic,'TooltipString','Draw square  Ctrl+2');
temppic = imread('arrow','bmp');
temppic = changepicbackground(temppic,bcolor);
toolbutton(3) = uicontrol('Style','togglebutton','Tag','arrowbutton','Position',[2+(buttonsep*2) figsize(4)-20 30 20],...
        'Callback','matclust(''arrowbutton_Callback'')','CData',temppic,'TooltipString','Selection tool  Ctrl+3');
temppic = imread('mag','bmp');
temppic = changepicbackground(temppic,bcolor);
toolbutton(4) = uicontrol('Style','togglebutton','Tag','magbutton','Position',[2+(buttonsep*3) figsize(4)-20 30 20],...
        'Callback','matclust(''magbutton_Callback'')','CData',temppic,'TooltipString','Magnify  Ctrl+4');
temppic = imread('hand','bmp');
temppic = changepicbackground(temppic,bcolor);
toolbutton(5) = uicontrol('Style','togglebutton','Tag','handbutton','Position',[2+(buttonsep*4) figsize(4)-20 30 20],...
        'Callback','matclust(''handbutton_Callback'')','CData',temppic,'TooltipString','Navigate  Ctrl+5');
 
temppic = imread('wandpic','bmp');
temppic = changepicbackground(temppic,bcolor);
toolbutton(6) = uicontrol('Style','togglebutton','Tag','wandbutton','Position',[2+(buttonsep*5) figsize(4)-20 30 20],...
        'Callback','matclust(''wandbutton_Callback'')','CData',temppic,'TooltipString','Magic wand  Ctrl+6','visible','on');
uicontrol('Style','slider','Tag','wandSlider','BackgroundColor',[.95 .95 .95],'ForegroundColor', figattrib.figureColor, ... 
    'Position',[7+(buttonsep*6) figsize(4)-22 200 20],'Max',100,'Min',10, ... 
    'SliderStep',[1 .1],'Value',25,'visible','off');
uicontrol('Style','text','Tag','wandSliderLabel','BackgroundColor',figattrib.figureColor, ... 
    'Position',[215+(buttonsep*6) figsize(4)-20 200 15],'string', 'Cluster size (few <--> many points)','HorizontalAlignment','left','visible','off');


set(toolbutton(1),'Value',1);
figattrib.toolselect = 1;

zoominbutton = uicontrol('Style','pushbutton','Tag','zoominbutton','Position',[figsize(3)-300 figsize(4)-20 30 20],...
        'Callback','matclust(''zoomIn'')','String','+','BackgroundColor',figattrib.figureColor);
zoomoutbutton = uicontrol('Style','pushbutton','Tag','zoomoutbutton','Position',[figsize(3)-330 figsize(4)-20 30 20],...
        'Callback','matclust(''zoomOut'')','String','-','BackgroundColor',figattrib.figureColor);
%*****************************************************
fprintf('.')
%create the right-side control panel (cluster controls, filter controls, etc)
%******************************************************
%create middle frame behind the cluster controls    
framecontrol(2) = uicontrol('Style','frame','Tag','frame2','BackgroundColor',figattrib.figureColor,'ForegroundColor', [0 0 0],'Position',[figsize(3)-260 1*figsize(4)/3 265 figsize(4)/3]);

%create the cluster controls
topyloc = (2*figsize(4)/3)-24;
for i = 1:length(figattrib.mixcolor)
    
    figattrib.clustcontrolcontext(i) = uicontextmenu('UserData',i);
    
    clustchoice(i,1) = uimenu(figattrib.clustcontrolcontext(i),'Tag',['clustcopyMenu',num2str(i)],'Label','Copy cluster to...','Callback','matclust(''copyclustMenuFcn'',get(gcbo,''UserData''),guidata(gcbo))','UserData',i,'Enable','off');
    clustchoice(i,2) = uimenu(figattrib.clustcontrolcontext(i),'Tag',['epochcopyMenu',num2str(i)],'Label','Copy all polygons in a time filter','Callback','matclust(''copyepochMenuFcn'',get(gcbo,''UserData''),guidata(gcbo))','UserData',i,'Enable','off');  
    clustchoice(i,3) = uimenu(figattrib.clustcontrolcontext(i),'Tag',['clustdeleteMenu',num2str(i)],'Label','Delete','UserData',i,'Enable','off');
    clustdeleteMenu(i,1) = uimenu(clustchoice(i,3),'Tag','clustdelallMenu','Label','Entire cluster','Callback','matclust(''clustdelete'',get(gcbo,''UserData''),guidata(gcbo),1)','UserData',i);
    clustdeleteMenu(i,2) = uimenu(clustchoice(i,3),'Tag','clustdelinfiltMenu','Label','Only in current time filter(s)','Callback','matclust(''clustfiltdelete'',get(gcbo,''UserData''),guidata(gcbo))','UserData',i);
    clustdeleteMenu(i,3) = uimenu(clustchoice(i,3),'Tag','clustdelwindowtMenu','Label','Polygon(s) in current window','Callback','matclust(''clustwindelete'',get(gcbo,''UserData''),guidata(gcbo))','UserData',i);
    clustdeleteMenu(i,4) = uimenu(clustchoice(i,3),'Tag','clustdelsingexcludesMenu','Label','Single point excludes','Callback','matclust(''spexcludedelete'',get(gcbo,''UserData''),guidata(gcbo))','UserData',i);
    clustchoice(i,4) = uimenu(figattrib.clustcontrolcontext(i),'Tag',['colorMenu',num2str(i)],'Label','Change color','Callback','getcolor(gcbo,guidata(gcbo))','UserData',i,'Enable','off');
    clustchoice(i,5) = uimenu(figattrib.clustcontrolcontext(i),'Tag',['orderMenu',num2str(i)],'Label','Order','UserData',i,'Enable','off');
    clustorderMenu(i,1) = uimenu(clustchoice(i,5),'Label','To Front','UserData',i,'Callback','matclust(''tofrontFcn'',gcbo,guidata(gcbo))');
    clustorderMenu(i,2) = uimenu(clustchoice(i,5),'Label','To Back','UserData',i,'Callback','matclust(''tobackFcn'',gcbo,guidata(gcbo))');
    clustorderMenu(i,3) = uimenu(clustchoice(i,5),'Label','Foreward One','UserData',i,'Callback','matclust(''forewardoneFcn'',gcbo,guidata(gcbo))');
    clustorderMenu(i,4) = uimenu(clustchoice(i,5),'Label','Back One','UserData',i,'Callback','matclust(''backoneFcn'',gcbo,guidata(gcbo))');
    clustchoice(i,6) = uimenu(figattrib.clustcontrolcontext(i),'Tag',['ClustToolsMenu',num2str(i)],'Label','Tools','UserData',i,'Enable','off');
    clusttoolschoice(i,1) = uimenu(clustchoice(i,6),'Tag',['viewISIMenu',num2str(i)],'Label','View ISI','Callback','ISIviewer(get(gcbo,''UserData''),guidata(gcbo))','UserData',i);
    clusttoolschoice(i,2) = uimenu(clustchoice(i,6),'Tag',['viewEvents',num2str(i)],'Label','View spikes','Callback','eventviewer(get(gcbo,''UserData''),guidata(gcbo))','UserData',i);
      
    % add any additional tool functions from the 'ClustTools' folder
	currdir = pwd;
	cd(figattrib.foldername);
	cd('ClustTools');
	toolsdir = dir;
	for j = 3:length(toolsdir)
        f = findstr(toolsdir(j).name,'.');
        tmpname = toolsdir(j).name(1:f-1);
        uimenu(clustchoice(i,6),'Label',tmpname,'Callback',['matclust(''ExecuteClustTool'',''',tmpname,''',get(gcbo,''UserData''))'],'UserData',i);
	end
	cd(currdir);
    if (mod(i,2))
        xloc = figsize(3)-242;
        yloc = topyloc-22*floor((i-1)/2)+1;
    else
        xloc = figsize(3)-180;
        yloc = topyloc-22*floor((i-1)/2)+1;    
    end
    
    figattrib.clustcontrol(i,1) = uicontrol('Style','frame','BackgroundColor',figattrib.figureColor,'ForegroundColor', [0 0 0],'Position',[xloc yloc-1 62 22],'UserData',i);
%    figattrib.clustcontrol(i,1) = uicontrol('Style','frame','BackgroundColor',figattrib.mixcolor(i,:),'ForegroundColor', [0 0 0],'Position',[xloc yloc-1 62 22],'UserData',i);
    
    figattrib.clustcontrol(i,2) = uicontrol('Style','togglebutton','BackgroundColor',figattrib.mixcolor(i,:),'Position',[xloc+1 yloc 20 20],'String',num2str(i), ...
        'Callback','matclust(''clustpickbutton_Callback'',gcbo,guidata(gcbo))','UserData',i,'UIContextMenu',figattrib.clustcontrolcontext(i),'TooltipString',['Cluster ',num2str(i)]);
    figattrib.clustcontrol(i,3) = uicontrol('Style','radiobutton','BackgroundColor',figattrib.mixcolor(i,:),'Position',[xloc+21 yloc 20 20],'Visible','off', ... 
        'UserData',i,'UIContextMenu',figattrib.clustcontrolcontext(i),'Callback','matclust(''clustradiobutton1_Callback'',gcbo,guidata(gcbo))','TooltipString','Filter out');
    figattrib.clustcontrol(i,4) = uicontrol('Style','radiobutton','BackgroundColor',figattrib.mixcolor(i,:),'Position',[xloc+41 yloc 20 20],'Visible','off', ... 
        'UserData',i,'UIContextMenu',figattrib.clustcontrolcontext(i),'Callback','matclust(''clustradiobutton2_Callback'',gcbo,guidata(gcbo))','Value',1,'TooltipString','Display');
end
set(figattrib.clustcontrol(1,2),'Value',1);
clustattrib.currclust = 1;
buttonlistlength = 11*length(figattrib.mixcolor);
uicontrol('Style','slider','Tag','clustslider','BackgroundColor',[.95 .95 .95],'ForegroundColor', figattrib.figureColor,'Position',[figsize(3)-259 figsize(4)/3 17 figsize(4)/3-1], ...
    'Max',buttonlistlength,'Min',20,'SliderStep',[(1/length(figattrib.mixcolor))*2 .1],'Value',buttonlistlength,'Callback','matclust(''clustslider_Callback'',gcbo,[],guidata(gcbo))');
figattrib.oldsliderval = buttonlistlength;

%create the cluster 0 button
clust0controlcontext = uicontextmenu;
clust0choice(1) = uimenu(clust0controlcontext,'Tag','color0Menu','Label','Change color','Callback','getcolor(gcbo,guidata(gcbo))','UserData',0);
clustattrib.cluster0button = uicontrol('Style','togglebutton','Tag','cluster0button','BackgroundColor',clustattrib.cluster0attrib.color,'Position',[figsize(3)-115 (2*figsize(4)/3)-23 110 20],'String','Cluster 0: visible', ...
        'Callback','matclust(''cluster0button_Callback'',gcbo,guidata(gcbo))','Value',1,'UIContextMenu',clust0controlcontext,'TooltipString','Cluster 0');
if (max(clustattrib.cluster0attrib.color) < .5)
    set(clustattrib.cluster0button,'ForegroundColor',[1 1 1]);
else
    set(clustattrib.cluster0button,'ForegroundColor',[0 0 0]);
end
fprintf('.')
%create the hide all and exclude all controls
hidealltext = uicontrol('Style','text','Tag','hidealltext','FontUnits','pixels','FontSize',11,'Position',[figsize(3)-118 2*(figsize(4)/3)-49 45 20],'String','Filter all','BackgroundColor',figattrib.figureColor);      
hideallonbutton = uicontrol('Style','pushbutton','Tag','hideallonbutton','FontUnits','pixels','FontSize',11,'Position',[figsize(3)-75 2*(figsize(4)/3)-45 25 17],'String','Set', ...
        'Callback','matclust(''hideallonbutton_Callback'',gcbo,guidata(gcbo))','BackgroundColor',figattrib.figureColor);
hidealloffbutton = uicontrol('Style','pushbutton','Tag','hidealloffbutton','FontUnits','pixels','FontSize',11,'Position',[figsize(3)-50 2*(figsize(4)/3)-45 45 17],'String','Release', ...
        'Callback','matclust(''hidealloffbutton_Callback'',gcbo,guidata(gcbo))','BackgroundColor',figattrib.figureColor);
excludealltext = uicontrol('Style','text','Tag','excludealltext','FontUnits','pixels','FontSize',11,'Position',[figsize(3)-118 2*(figsize(4)/3)-69 46 20],'String','Display','BackgroundColor',figattrib.figureColor);      
excludeallonbutton = uicontrol('Style','pushbutton','Tag','excludeallonbutton','FontUnits','pixels','FontSize',11,'Position',[figsize(3)-75 2*(figsize(4)/3)-65 25 17],'String','Set', ...
        'Callback','matclust(''excludeallonbutton_Callback'',gcbo,guidata(gcbo))','BackgroundColor',figattrib.figureColor);
excludealloffbutton = uicontrol('Style','pushbutton','Tag','excludealloffbutton','FontUnits','pixels','FontSize',11,'Position',[figsize(3)-50 2*(figsize(4)/3)-65 45 17],'String','Release', ...
        'Callback','matclust(''excludealloffbutton_Callback'',gcbo,guidata(gcbo))','BackgroundColor',figattrib.figureColor);
    
clusterinfolist = uicontrol('Style','listbox','Tag','clusterinfo','Position',[figsize(3)-115 (figsize(4)/3)+3 110 (figsize(4)/3)-75], ...
        'BackgroundColor',[1 1 1],'Max',100,'Min',0,'Callback','matclust(''clusterinfo_Callback'',gcbo,guidata(gcbo))');
    
    
%create the top and bottom frames in the controls panel
framecontrol(1) = uicontrol('Style','frame','Tag','frame1','BackgroundColor',figattrib.figureColor,'ForegroundColor', [0 0 0],'Position',[figsize(3)-260 (2*figsize(4)/3)-1 265 (figsize(4)/3)+10]);
framecontrol(3) = uicontrol('Style','frame','Tag','frame3','BackgroundColor',figattrib.figureColor,'ForegroundColor', [0 0 0],'Position',[figsize(3)-260 20 265 (figsize(4)/3)-19]);

%create the listboxes to control the current clustering axes
axeslistcontext = uicontextmenu('Tag','axeslistContext1');
acontext1 = uimenu(axeslistcontext,'Tag','paramaddContext1','Label','Add rotation','Callback','matclust(''paramrotate'',gcbo,guidata(gcbo))');
acontext2 = uimenu(axeslistcontext,'Tag','paramdeleteContext1','Label','Delete rotation','Callback','matclust(''paramrotationdelete'',gcbo,guidata(gcbo))','UserData',1);
acontext3 = uimenu(axeslistcontext,'Tag','parameditContext1','Label','Edit name','Callback','matclust(''parameditname'',gcbo,guidata(gcbo))','UserData',1);

set(acontext1,'Enable','off');
set(acontext2,'Enable','off');
set(acontext3,'Enable','off');

listcontrol(1) = uicontrol('Style','listbox','Tag','listbox1','FontUnits','pixels','FontSize',11,'BackgroundColor',[1 1 1],'ForegroundColor',[0 0 0], ...
'Position',[figsize(3)-255 (2*figsize(4)/3)+1 120 (figsize(4)/3)-20],'String',clustdata.names,'Max',1,'Min',1, ...
'Callback','matclust(''listbox1_Callback'',gcbo,[],guidata(gcbo))', 'Value',[1],'UIContextMenu',axeslistcontext,'Enable','off');

axeslistcontext = uicontextmenu('Tag','axeslistContext2');
acontext1 = uimenu(axeslistcontext,'Tag','paramaddContext2','Label','Add rotation','Callback','matclust(''paramrotate'',gcbo,guidata(gcbo))');
acontext2 = uimenu(axeslistcontext,'Tag','paramdeleteContext2','Label','Delete rotation','Callback','matclust(''paramrotationdelete'',gcbo,guidata(gcbo))','UserData',2);
acontext3 = uimenu(axeslistcontext,'Tag','parameditContext2','Label','Edit name','Callback','matclust(''parameditname'',gcbo,guidata(gcbo))','UserData',2);

set(acontext1,'Enable','off');
set(acontext2,'Enable','off');
set(acontext3,'Enable','off');


listcontrol(2) = uicontrol('Style','listbox','Tag','listbox2','FontUnits','pixels','FontSize',11,'BackgroundColor',[1 1 1],'ForegroundColor',[0 0 0], ...
'Position',[figsize(3)-125 (2*figsize(4)/3)+1 120 (figsize(4)/3)-20],'String',clustdata.names,'Max',1,'Min',1, ...
'Callback','matclust(''listbox2_Callback'',gcbo,[],guidata(gcbo))', 'Value',1,'UIContextMenu',axeslistcontext,'Enable','off');
set(listcontrol(1),'UserData',1);
set(listcontrol(2),'UserData',2);
xtext = uicontrol('Style','text','Tag','xText','Position',[figsize(3)-255 (figsize(4))-19 120 15],'String','X Axis','BackgroundColor',figattrib.figureColor);
ytext = uicontrol('Style','text','Tag','yText','Position',[figsize(3)-125 (figsize(4))-19 120 15],'String','Y Axis','BackgroundColor',figattrib.figureColor);

%create the filter listboxes
timefiltertext = uicontrol('Style','text','Tag','timefilterText','FontUnits','pixels','FontSize',11,'Position',[figsize(3)-255 (figsize(4)/3)-22 120 20],'String','Time','BackgroundColor',figattrib.figureColor);      
otherfiltertext = uicontrol('Style','text','Tag','otherfilterText','FontUnits','pixels','FontSize',11,'Position',[figsize(3)-125 (figsize(4)/3)-22 120 20],'String','Filters','BackgroundColor',figattrib.figureColor);      

timeFiltOptionButton = uicontrol('Style','radiobutton','Tag','timeFiltOptButton','Position',[figsize(3)-255 (figsize(4)/3)-16 32 17], ...
        'Callback','matclust(''TFiltBehaveFcn'',gcbo,guidata(gcbo))','Value',1,'UserData',3,'TooltipString','Time selection mode: one at a time','BackgroundColor',figattrib.figureColor);


%figattrib.clustcontrolcontext(i) = uicontextmenu('UserData',i);
%clustchoice(i,1) = uimenu(figattrib.clustcontrolcontext(i),'Tag',['colorMenu',num2str(i)],'Label','Change color','Callback','getcolor(gcbo,guidata(gcbo))','UserData',i);
timefiltercontext = uicontextmenu('Tag','timefilterContext');
tcontext1 = uimenu(timefiltercontext,'Tag','timeaddContext','Label','Add time filter','Callback','timefilteradd(gcbo,guidata(gcbo))');
tcontext2 = uimenu(timefiltercontext,'Tag','timedeleteContext','Label','Delete time filter','Callback','matclust(''timefilterdelete'',gcbo,guidata(gcbo))');
tcontext3 = uimenu(timefiltercontext,'Tag','timeeditContext','Label','Edit time filter','Callback','timefilteredit(gcbo,guidata(gcbo))');

set(tcontext1,'Enable','off');
set(tcontext2,'Enable','off');

otherfiltercontext = uicontextmenu('Tag','otherfilterContext');
ocontext1 = uimenu(otherfiltercontext,'Tag','otheraddContext','Label','Add filter','Callback','otherfilteradd(gcbo,guidata(gcbo))');
ocontext2 = uimenu(otherfiltercontext,'Tag','otherdeleteContext','Label','Delete filter','Callback','matclust(''otherfilterdelete'',gcbo,guidata(gcbo))');
set(ocontext1,'Enable','off');
set(ocontext2,'Enable','off');
for bfill = 1:32
    blanks(bfill) = {[num2str(bfill)]};
end
blanks{1} = [num2str(1),187,' All points'];

clustdata.timefilternames = blanks;
clustdata.otherfilternames = blanks;
clustdata.timefilterranges = clustdata.datarange(:,1)';
listcontrol(3) = uicontrol('Style','listbox','Tag','timefilterList','FontUnits','pixels','FontSize',11,'BackgroundColor',[1 1 1],'ForegroundColor',[0 0 0], ...
    'Position',[figsize(3)-255 25 120 (figsize(4)/3)-40],'String',blanks,'Max',2,'Min',0,'Value',[], ...
    'Callback','matclust(''timefilterList_Callback'',gcbo,guidata(gcbo))','UIContextMenu',timefiltercontext,'UserData',[1 0],'Enable','off');
listcontrol(4) = uicontrol('Style','listbox','Tag','otherfilterList','FontUnits','pixels','FontSize',11,'BackgroundColor',[1 1 1],'ForegroundColor',[0 0 0], ...
    'Position',[figsize(3)-125 25 120 (figsize(4)/3)-40],'String',blanks,'Max',2,'Min',0,'Value',[], ...
    'Callback','matclust(''otherfilterList_Callback'',gcbo,guidata(gcbo))','UIContextMenu',otherfiltercontext,'UserData',[1 0],'Enable','off');
%**************************************************




%create the status box
statusbox = uicontrol('Style','frame','Tag','statusBox','BackgroundColor',figattrib.figureColor,'ForegroundColor', [1 1 1],'Position',[0 0 figsize(3)+5 20]);
statustext = uicontrol('Style','text','Tag','statusText','FontUnits','pixels','FontSize',11,'Position',[0 0 figsize(3) 18],'String','  Ready','BackgroundColor',figattrib.figureColor,'HorizontalAlignment','left','UserData','Ready');      

fprintf('\n')
set(h,'Visible','on')

%create timer controls
figattrib.timers = timer('Period',.5,'TimerFcn','matclust(''redrawTimerFcn'')','TasksToExecute',100,'ExecutionMode','fixedDelay');

% Generate a structure of handles to pass to callbacks, and store it. 
handles = guihandles(h);
handles.clustcontrol = figattrib.clustcontrol;
handles.polycontext = polycontext;
handles.mainClusterMenu = mainClusterMenu;
figattrib.handles = handles;
guidata(h, handles);

%execute sequence of startup commands 
%*********************************
plotgraph(handles);
currdir = pwd;
cd(figattrib.datafoldername);
mkdir M_dataopen
cd(figattrib.datafoldername);
cd('M_dataopen');
dirnames = dir;

if (length(dirnames)>2)
       
	cd(figattrib.datafoldername);
	if (ispc)
        system('rd /S /Q M_dataopen');
	elseif (isunix)
        system('rm -r -f M_dataopen');
	end
	mkdir('M_dataopen');         
end
cd(currdir);
%*********************************

%---------------------------------------------------------------------




%Section for callbacks of buttons and lists, etc
%---------------------------------------------------------------------
% --------------------------------------------------------------------
function varargout = listbox1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.listbox1.

global graphattrib;
global clustdata;
global figattrib;

v = (get(h,'Value'));
if (clustdata.filledparam(v))
    releasepolydraw(handles);
	hidecurrentpolygons(handles);
    figattrib.cashedPolygon = [];
    set(handles.pastePolygonFilterMainmenu,'Enable','off');
	
	%change the axes
	a1 = get(handles.listbox1,'Value');
	a2 = get(handles.listbox2,'Value');
	set(handles.listbox1,'UserData',a1);
	set(handles.listbox2,'UserData',a2);
	
	%show the new polygons
	showcurrentpolygons(handles);
    
	if (a1>clustdata.origparams)
        set(handles.paramdeleteContext1,'Enable','on');
        set(handles.paramdeleteContext2,'Enable','on');
        a2 = a1;
        set(handles.listbox2,'UserData',a2);
        set(handles.listbox2,'Value',a2);
        graphattrib.currentViewIsRotation = 1;
	else
        set(handles.paramdeleteContext1,'Enable','off');
        set(handles.paramdeleteContext2,'Enable','off');
        try
            CVR = graphattrib.currentViewIsRotation;
        catch
            CVR = 0;
        end
        if CVR
           a2 = a1;
           set(handles.listbox2,'UserData',a2);
           set(handles.listbox2,'Value',a2);      
        end
        graphattrib.currentViewIsRotation = 0;
	end
	updateclustinfo(handles);
	plotgraph(handles)
else
    a1 = get(handles.listbox1,'UserData');
	a2 = get(handles.listbox2,'UserData');
    set(handles.listbox1,'Value',a1);
    set(handles.listbox2,'Value',a2);
end

% --------------------------------------------------------------------
function varargout = listbox2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.listbox2.

global graphattrib;
global clustdata;
global figattrib;

v = (get(h,'Value'));
if (clustdata.filledparam(v))
	releasepolydraw(handles);
    hidecurrentpolygons(handles);
    figattrib.cashedPolygon = [];
    set(handles.pastePolygonFilterMainmenu,'Enable','off');
	
	%change the axes
	a1 = get(handles.listbox1,'Value');
	a2 = get(handles.listbox2,'Value');
	set(handles.listbox1,'UserData',a1);
	set(handles.listbox2,'UserData',a2);
	
	%show the new polygons
	showcurrentpolygons(handles);
    
	if (a2>clustdata.origparams)
        set(handles.paramdeleteContext1,'Enable','on');
        set(handles.paramdeleteContext2,'Enable','on');
        a1 = a2;
        set(handles.listbox1,'UserData',a1);
        set(handles.listbox1,'Value',a1);
        graphattrib.currentViewIsRotation = 1;
        
        
    else
        set(handles.paramdeleteContext1,'Enable','off');
        set(handles.paramdeleteContext2,'Enable','off');
        try
            CVR = graphattrib.currentViewIsRotation;
        catch
            CVR = 0;
        end
        if CVR
           a1 = a2;
           set(handles.listbox1,'UserData',a1);
           set(handles.listbox1,'Value',a1);
        end
        graphattrib.currentViewIsRotation = 0;
        
	end
	updateclustinfo(handles);
	plotgraph(handles)
else
    a1 = get(handles.listbox1,'UserData');
	a2 = get(handles.listbox2,'UserData');
    set(handles.listbox1,'Value',a1);
    set(handles.listbox2,'Value',a2);
end

% ------------------------------------------------------------
function polygonbutton_Callback()

global figattrib;
global graphattrib;
handles = figattrib.handles;
hObject = handles.polygonbutton;
clearhighlight;
figattrib.toolselect = 1;
set(handles.figure1,'Pointer','crosshair');
if (get(hObject,'Value')==0)
    set(hObject,'Value',1);
end

resetButtonBackgrounds(handles)
bcolor = [200 200 200];
temppic = imread('poly','bmp');
temppic = changepicbackground(temppic,bcolor);
set(handles.polygonbutton,'CData',temppic);
set(handles.squarebutton,'Value',0);
set(handles.arrowbutton,'Value',0);
set(handles.magbutton,'Value',0);
set(handles.handbutton,'Value',0);
set(handles.wandbutton,'Value',0);
set(handles.wandSlider,'visible','off');
set(handles.wandSliderLabel,'visible','off');
% -------------------------------------------------------------
function squarebutton_Callback()

global figattrib;
global graphattrib;
handles = figattrib.handles;

hObject = handles.squarebutton;
releasepolydraw(handles);
clearhighlight;
figattrib.toolselect = 2;
set(handles.figure1,'Pointer','crosshair');
if (get(hObject,'Value')==0)
    set(hObject,'Value',1);
end

resetButtonBackgrounds(handles)
bcolor = [200 200 200];
temppic = imread('square','bmp');
temppic = changepicbackground(temppic,bcolor);
set(handles.squarebutton,'CData',temppic);

set(handles.polygonbutton,'Value',0);
set(handles.arrowbutton,'Value',0);
set(handles.magbutton,'Value',0);
set(handles.handbutton,'Value',0);
set(handles.wandbutton,'Value',0);
set(handles.wandSlider,'visible','off');
set(handles.wandSliderLabel,'visible','off');
% --------------------------------------------------------------
function arrowbutton_Callback()

global figattrib;
global graphattrib;
handles = figattrib.handles;

hObject = handles.arrowbutton;
releasepolydraw(handles);
figattrib.toolselect = 3;
set(handles.figure1,'Pointer','arrow');
if (get(hObject,'Value')==0)
    set(hObject,'Value',1);
end

resetButtonBackgrounds(handles)
bcolor = [200 200 200];
temppic = imread('arrow','bmp');
temppic = changepicbackground(temppic,bcolor);
set(handles.arrowbutton,'CData',temppic);
set(handles.polygonbutton,'Value',0);
set(handles.squarebutton,'Value',0);
set(handles.magbutton,'Value',0);
set(handles.handbutton,'Value',0);
set(handles.wandbutton,'Value',0);
set(handles.wandSlider,'visible','off');
set(handles.wandSliderLabel,'visible','off');
% --------------------------------------------------------------
function magbutton_Callback()

global figattrib;
global graphattrib;
handles = figattrib.handles;

hObject = handles.magbutton;

releasepolydraw(handles);
clearhighlight;
figattrib.toolselect = 4;
load magpic;
set(handles.figure1,'Pointer','custom','PointerShapeCData',magpic,...
            'PointerShapeHotSpot',[6 5]);
if (get(hObject,'Value')==0)
    set(hObject,'Value',1);
end

resetButtonBackgrounds(handles)
bcolor = [200 200 200];
temppic = imread('mag','bmp');
temppic = changepicbackground(temppic,bcolor);
set(handles.magbutton,'CData',temppic);
set(handles.polygonbutton,'Value',0);
set(handles.squarebutton,'Value',0);
set(handles.arrowbutton,'Value',0);
set(handles.handbutton,'Value',0);
set(handles.wandbutton,'Value',0);
set(handles.wandSlider,'visible','off');
set(handles.wandSliderLabel,'visible','off');
%------------------------------------------------------------
function handbutton_Callback()

global figattrib;
global graphattrib;
handles = figattrib.handles;

hObject = handles.handbutton;
releasepolydraw(handles);
clearhighlight;
figattrib.toolselect = 5;
load handpointer;
set(handles.figure1,'Pointer','custom','PointerShapeCData',handpointer,...
            'PointerShapeHotSpot',[8 8]);
if (get(hObject,'Value')==0)
    set(hObject,'Value',1);
end

resetButtonBackgrounds(handles)
bcolor = [200 200 200];
temppic = imread('hand','bmp');
temppic = changepicbackground(temppic,bcolor);
set(handles.handbutton,'CData',temppic);
set(handles.polygonbutton,'Value',0);
set(handles.squarebutton,'Value',0);
set(handles.arrowbutton,'Value',0);
set(handles.magbutton,'Value',0);
set(handles.wandbutton,'Value',0);
set(handles.wandSlider,'visible','off');
set(handles.wandSliderLabel,'visible','off');
%----------------------------------------------------------------
function wandbutton_Callback()

global figattrib;
global graphattrib;
handles = figattrib.handles;

hObject = handles.wandbutton;
clearhighlight;
figattrib.toolselect = 6;
load  wandpic;
set(handles.figure1,'Pointer','custom','PointerShapeCData',wandpic,...
            'PointerShapeHotSpot',[6 6]);
if (get(hObject,'Value')==0)
    set(hObject,'Value',1);
end

resetButtonBackgrounds(handles)
bcolor = [200 200 200];
temppic = imread('wandpic','bmp');
temppic = changepicbackground(temppic,bcolor);
set(handles.wandbutton,'CData',temppic);
set(handles.polygonbutton,'Value',0);
set(handles.squarebutton,'Value',0);
set(handles.arrowbutton,'Value',0);
set(handles.magbutton,'Value',0);
set(handles.handbutton,'Value',0);
set(handles.wandSlider,'visible','on');
set(handles.wandSliderLabel,'visible','on');
%----------------------------------------------

function clustslider_Callback(hObject, eventdata, handles)
%controls the slider and visible cluster controls

global figattrib;

newval = get(hObject,'Value');
change = figattrib.oldsliderval - newval;
currloc = get(figattrib.clustcontrol,'Position');
for i = 1:length(currloc)
    currloc{i}(2) = currloc{i}(2)+change;
end
set(figattrib.clustcontrol,{'Position'},currloc);
figattrib.oldsliderval = newval;
%-------------------------------------------------------------
function clustpickbutton_Callback(hObject,handles)
%called when any of the cluster number buttons are pressed

global figattrib;
global clustattrib;

releasepolydraw(handles);
clearhighlight;

set(figattrib.clustcontrol(clustattrib.currclust,2),'BackgroundColor',figattrib.mixcolor(clustattrib.currclust,:));
clustattrib.currclust = get(hObject,'UserData');
if (ismember(clustattrib.currclust,clustattrib.clustersOn))
    set(handles.mainClustMenu,'Enable','on');
    for menuNum = 1:length(handles.mainClusterMenu)
        set(handles.mainClusterMenu(menuNum),'UserData',clustattrib.currclust)
    end
else
    set(handles.mainClustMenu,'Enable','off');
end

updateclustinfo(handles);

set(figattrib.clustcontrol(:,2),'Value',0);
set(hObject,'Value',1);
    
try %if a dataset is loaded, then display depend info
    set(figattrib.handles.statusText,'String', ...
    ['Ready                                       CLUSTER ',num2str(clustattrib.currclust), ... 
    '-  Filtered clusters: ',num2str(find(clustattrib.dependencies(clustattrib.currclust,:))), ... 
    '    Dependent clusters: ', num2str(find(clustattrib.dependencies(:,clustattrib.currclust))'),'']);
    set(figattrib.handles.statusText,'UserData', 'Ready');
end
%------------------------------------------------------------
function clustradiobutton1_Callback(hObject,handles)
%the 'hide' radiobutton on the cluster controls

global clustattrib;

clustfocus = get(hObject,'UserData');
clustattrib.hiddenclusters(clustfocus,1) = get(hObject,'Value');
plotgraph(handles)
%------------------------------------------------------------
function clustradiobutton2_Callback(hObject,handles)
%the  'exclude' radiobutton on the cluster controls

global clustattrib;

clustfocus = get(hObject,'UserData');
clustattrib.hiddenclusters(clustfocus,2) = 1-get(hObject,'Value');
plotgraph(handles)
%-------------------------------------------------------------
function cluster0button_Callback(hObject,handles)
%toggles the visibility of cluster 0 when the 'cluster 0' button is pressed

global clustattrib;

clustattrib.cluster0attrib.show = get(hObject,'Value');
if (get(hObject,'Value'))
    set(hObject,'String','Cluster 0: visible');
else
    set(hObject,'String','Cluster 0: hidden');
end
plotgraph(handles)
%-------------------------------------------------------------
function hideallonbutton_Callback(hObject,handles)

global clustattrib;
global figattrib;

hidevalue = 1;
for i = 1:size(figattrib.clustcontrol,1)
     clusthideval(i,1) = {hidevalue};
end
set(figattrib.clustcontrol(:,3),{'Value'},clusthideval);
clustattrib.hiddenclusters(:,1) = 1;
plotgraph(handles)
%-------------------------------------------------------------
function hidealloffbutton_Callback(hObject,handles)

global clustattrib;
global figattrib;

hidevalue = 0;
for i = 1:size(figattrib.clustcontrol,1)
     clusthideval(i,1) = {hidevalue};
end
set(figattrib.clustcontrol(:,3),{'Value'},clusthideval);
clustattrib.hiddenclusters(:,1) = 0;
plotgraph(handles)
%-------------------------------------------------------------
function excludeallonbutton_Callback(hObject,handles)

global clustattrib;
global figattrib;

excludevalue = 1;
for i = 1:size(figattrib.clustcontrol,1)
     clustexcludeval(i,1) = {excludevalue};
end
set(figattrib.clustcontrol(:,4),{'Value'},clustexcludeval);
clustattrib.hiddenclusters(:,2) = 0;
plotgraph(handles)
%-------------------------------------------------------------
function excludealloffbutton_Callback(hObject,handles)

global clustattrib;
global figattrib;

excludevalue = 0;
for i = 1:size(figattrib.clustcontrol,1)
     clustexcludeval(i,1) = {excludevalue};
end
set(figattrib.clustcontrol(:,4),{'Value'},clustexcludeval);
clustattrib.hiddenclusters(:,2) = 1;
plotgraph(handles)
%-------------------------------------------------------------
function timefilterList_Callback(hObject,handles)

global clustdata;
global graphattrib;
global figattrib;

allowMultipleTimes = figattrib.tFiltAllowMultiple;
newval = get(hObject,'Value');
if ~isempty(newval)
	if isempty(find(clustdata.timefiltermemmap == newval))
        set(handles.timeaddContext,'Enable','on');
        set(handles.timedeleteContext,'Enable','off');
        set(handles.timeeditContext,'Enable','off');
	else
        set(handles.timeaddContext,'Enable','off');
        set(handles.timedeleteContext,'Enable','on');
        set(handles.timeeditContext,'Enable','on');
	end
	if (newval == 1)
        set(handles.timeaddContext,'Enable','off');
        set(handles.timedeleteContext,'Enable','off');
        set(handles.timeeditContext,'Enable','off');
	end
	
	names = get(hObject,'String');
	d = get(hObject,'UserData');
	set(hObject,'Userdata',[newval toc]);
	
    %if ((d(1) == newval)&(abs(toc-d(2))<.5)& ~isempty(find(clustdata.timefiltermemmap == newval))) %item was double-clicked
	if ((~isempty(find(clustdata.timefiltermemmap == newval))) && ~figattrib.tFiltAllowMultiple) || ... %single-click highlight when no multiples allowed 
            ((d(1) == newval)&(abs(toc-d(2))<.5)& ~isempty(find(clustdata.timefiltermemmap == newval))) %item was double-clicked
        releasepolydraw(handles);
        a1 = get(handles.listbox1,'UserData');
        a2 = get(handles.listbox2,'UserData');
        
        hidecurrentpolygons(handles); %hide the currently displayed polygons
        
        %if the user selected all times, then all other time filters are
        %de-activiated
        if (newval == 1)
            allowMultipleTimes = 0;
        end
        if (allowMultipleTimes)  %if multiple time filters are allowed on at the same time
            arrowcheck = strfind(names{newval},187);
            if ((~isempty(arrowcheck))) %turn off filter
                names{newval}(arrowcheck) = ' ';
                clustdata.timefiltersOn(find(clustdata.timefiltermemmap == newval)) = 0;    
            else %turn on filter
                spacecheck = min(strfind(names{newval},' '));
                names{newval}(spacecheck) = 187; %add double arrows to highlight selection
                clustdata.timefiltersOn(find(clustdata.timefiltermemmap == newval)) = 1;        
                names{1} = [num2str(1),'  All points'];
                clustdata.timefiltersOn(1) = 0;
            end
            %if all time filters are turned off, then turn on all points
            %filter
            if (length(find(clustdata.timefiltersOn)) == 0)
                names{1} = [num2str(1),187,' All points'];
                clustdata.timefiltersOn(1) = 1;
            end
                
        else %multiple time filters are not allowed
            clustdata.timefiltersOn = zeros(32,1)';
            clustdata.timefiltersOn(find(clustdata.timefiltermemmap == newval)) = 1; %finds the index to the filter information
            for i = 1:length(names)
                tempname = names{i};
                tempname(strfind(names{i},187)) = ' ';
                
                if (i == newval)
                    tempname(min(strfind(tempname,' '))) = 187; %add double arrows to highlight selection
                end
                names{i} = tempname;
            end
            
        end
        
        
        set(hObject,'String',names);
        clustdata.timefilternames = names;

        FilterPoints;
        	
        showcurrentpolygons(handles); %show the new polygons
        
        %cleartimefilter;
        updateclustinfo(handles);
        plotgraph(handles)
        %addnewstate('filter points',handles); 
	end
end
%--------------------------------------------------------------------------
function otherfilterList_Callback(hObject,handles)

global clustdata;
global graphattrib;
global figattrib;

clearhighlight;
newval = get(hObject,'Value');
if ~isempty(newval)
	if isempty(find(clustdata.otherfiltermemmap == newval))
        set(handles.otheraddContext,'Enable','on');
        set(handles.otherdeleteContext,'Enable','off');
	else
        set(handles.otheraddContext,'Enable','off');
        set(handles.otherdeleteContext,'Enable','on');
	end
	if (newval == 1)
        set(handles.otheraddContext,'Enable','off');
        set(handles.otherdeleteContext,'Enable','off');
	end
	
	names = get(hObject,'String');
	d = get(hObject,'UserData');
	set(hObject,'Userdata',[newval toc]);
	if ((d(1) == newval)&(abs(toc-d(2))<.5)& ~isempty(find(clustdata.otherfiltermemmap == newval)) & (newval~= 1)) %item was double-clicked
		releasepolydraw(handles);
        
        %figattrib.cashedPolygon = [];
        %set(handles.pastePolygonFilterMainmenu,'Enable','off');
        a1 = get(handles.listbox1,'UserData');
        a2 = get(handles.listbox2,'UserData');
        
        hidecurrentpolygons(handles); %hide the currently displayed polygons
        
        arrowcheck = strfind(names{newval},187);
        if ~isempty(arrowcheck)
            names{newval}(arrowcheck) = ' ';
            clustdata.otherfiltersOn(find(clustdata.otherfiltermemmap == newval)) = 0;    
        else
            spacecheck = min(strfind(names{newval},' '));
            names{newval}(spacecheck) = 187; %add double arrows to highlight selection
            clustdata.otherfiltersOn(find(clustdata.otherfiltermemmap == newval)) = 1;        
        end
        
        set(hObject,'String',names);
        
        clustdata.otherfilternames = names;
        FilterPoints;
        showcurrentpolygons(handles); %show the new polygons
        
        clearotherfilter;
        updateclustinfo(handles);
        plotgraph(handles)
        %addnewstate('filter points',handles); 
	end
end
%--------------------------------------------------------------------------
function clusterinfo_Callback(hObject,handles)

global clustattrib;
global clustdata;
global graphattrib;


releasepolydraw(handles);
hidecurrentpolygons(handles);

if (clustattrib.nodata)
    return  %if no data is loaded, then do nothing
end

currval = get(hObject,'Value');
axes = clustattrib.clusters{clustattrib.currclust}.defineaxes(currval,:);
a1 = axes(1);
a2 = axes(2);
filters = clustattrib.filterindex{axes(3)};
t = filters(1);
o = filters(2:end);

names = get(handles.timefilterList,'String');

clustdata.timefiltersOn = zeros(32,1)';

for i = 1:length(names)
    tempname = names{i};
    tempname(strfind(names{i},187)) = ' ';
       
    if ismember(i,t)
        tempname(min(strfind(tempname,' '))) = 187; %add double arrows to highlight selection
        clustdata.timefiltersOn(i) = 1; %[clustdata.otherfiltersOn find(clustdata.otherfiltermemmap == i)];   
    end
    names{i} = tempname;
end

set(handles.timefilterList,'String',names);

names = get(handles.otherfilterList,'String');
clustdata.otherfiltersOn = zeros(32,1)';
clustdata.otherfiltersOn(1) = 1;
for i = 1:length(names)
    tempname = names{i};
    tempname(strfind(names{i},187)) = ' ';
       
    if ismember(i,o)
        tempname(min(strfind(tempname,' '))) = 187; %add double arrows to highlight selection
        clustdata.otherfiltersOn(i) = 1; %[clustdata.otherfiltersOn find(clustdata.otherfiltermemmap == i)];   
    end
    names{i} = tempname;
end
set(handles.otherfilterList,'String',names);


FilterPoints;
set(handles.listbox1,'UserData',a1);
set(handles.listbox2,'UserData',a2);
set(handles.listbox1,'Value',a1);
set(handles.listbox2,'Value',a2);

showcurrentpolygons(handles);   
%cleartimefilter;
updateclustinfo(handles);
plotgraph(handles)
%----------------------------------------------
function deletepolymenuCallback(hObject, handles)

global graphattrib;

%info = get(graphattrib.currentpolyhighlight,'UserData');
%deletepoly(info(1),info(2),info(3),info(4));
%updateclustinfo(handles);    
%plotgraph(handles);
%addnewstate('delete polygon',handles); 

polyhandle = graphattrib.currentpolyhighlight;
totalpolyindex = get(polyhandle,'UserData');
if ~iscell(totalpolyindex)
    tmppolyindex{1} = totalpolyindex;
    totalpolyindex = tmppolyindex;
end
for i = 1:length(totalpolyindex)
    polyindex = totalpolyindex{i};
    deletepoly(polyindex(1),polyindex(2),polyindex(3),polyindex(4));
end
updateclustinfo(handles);
plotgraph(handles);
addnewstate('delete polygon',handles);
%---------------------------------------------
function addPolygonFilter_Callback(hObject, handles)

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;

currentPoly = get(graphattrib.currentpolyhighlight,'UserData');
a1 = currentPoly(1);
a2 = currentPoly(2);
findex = currentPoly(3);
clustnum = currentPoly(4);
vertices = graphattrib.polyg(a1,a2,findex).relativevertices{clustnum}; %[x y] points of the polygon (n by 2 matrix)

[filename, path] = filebrowse('open','filter','*.m','title','Choose filter function');
if (isempty(filename))
    return
end
funcname = filename(1:strfind(filename,'.')-1);

filts = clustattrib.filterindex{findex};
tmpfilter = fastandbit(clustdata.timefilters,filts(1)); %temporal filters
tmpfilter2 = fastandbit(clustdata.otherfilters,filts(2:end)); %other filters
filteredpoints = tmpfilter & tmpfilter2;
in = find(filteredpoints & graphattrib.nonhiddenpoints); %don't consider currently hidden clusters
senddata = clustdata.params(in,:); %only the currently filtered data is sent to the program

%scale the polygon to the data range so that the polygon matches up
currentdatarange = [clustdata.datarange(1,a1);clustdata.datarange(2,a1)];
vertices(:,1) = (vertices(:,1)*(currentdatarange(2)-currentdatarange(1)))+currentdatarange(1);
currentdatarange = [clustdata.datarange(1,a2);clustdata.datarange(2,a2)];
vertices(:,2) = (vertices(:,2)*(currentdatarange(2)-currentdatarange(1)))+currentdatarange(1);


%scale the data to the range (so that the polygon matches up
% for i = 1:size(clustdata.params,2)    
%     currentdatarange = [clustdata.datarange(1,i);clustdata.datarange(2,i)];    
%     senddata(:,i) = (senddata(:,i) - currentdatarange(1))/(currentdatarange(2)-currentdatarange(1));
% end

currdir = pwd;
cd(path);
eval(['funchandle = @',funcname,';']);
cd(currdir);

result = [];
%The function has three inputs 1)the polygon x and y vertices, 2) the 2 dimension
%indices of the polygon's x and y vertices, 3) the full multidimensional data (n by numDims).
%The output is a vector of 0's and 1's the same length as n.

result = feval(funchandle,vertices, [a1 a2], senddata); 

result = (result(:) ~= 0);
if (length(result) ~= size(senddata,1))
    errordlg('The output of the filter function must be a vector of the same length as the input data.','Error');
    return
end

polynum = clustattrib.clusters{clustnum}.polyindex( ...
   find((clustattrib.clusters{clustnum}.polyindex(:,1)==a1)&(clustattrib.clusters{clustnum}.polyindex(:,2)==a2)&(clustattrib.clusters{clustnum}.polyindex(:,3)==findex)),4);
memvector = ceil(polynum/32); 

excluding = false(length(clustdata.params(:,1)),1);
including = false(length(clustdata.params(:,1)),1); %default assumption: all points are outside the polygon
excluding(in) = true;  %default assumption: all points in the current filter are excluded by this polygon (note: points not in the filter are not excluded)
in = in(find(result)); %find the points of the current filter that are inside the polygon's filter

excluding(in) = false;  %points within polygon are not excluded
including(in) = true; %points within polygon are included

paramlength = length(clustdata.params(:,1));

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
clustattrib.clusters{clustnum}.index = uint32(find(~in & in2 & ~in3 & graphattrib.nonhiddenpoints));

%any clusters that are dependent on this cluster must also be redefined
%(this will work its way down the tree of dependencies)
for i = find(clustattrib.dependencies(:,clustnum))'
    clustexcludes = find(clustattrib.dependencies(i,:));
    findInPoints([],[],[],i, clustexcludes);
end

set(graphattrib.polyg(a1,a2,findex).lines{clustnum},'LineStyle','-.','LineWidth',2); 
graphattrib.polyg(a1,a2,findex).type{clustnum} = 2; %polygon is type 2, meaning it has a filter attached

plotgraph(handles);  
addnewstate('create polygon filter',handles); 
%---------------------------------------------

%Section for functions controlling the graph and objects within the graph
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function figure1_WindowButtonDownFcn(hObject, handles)
%called whenever the left mouse button in pressed within the figure 

global figattrib;
global graphattrib;
global clustattrib;
global clustdata;

a1 = get(handles.listbox1,'value');
a2 = get(handles.listbox2,'value');
[findex, filterline] = findfilterindex(handles);
findex = findex(1);
filterline = filterline(1,:);
[findex2, filterline2] = findfilterindex(handles);
stopwatch = toc;
tic;

figloc = get(graphattrib.graphwindow,'Position');
windowsize = get(handles.figure1,'Position');
M = get(hObject,'CurrentPoint');

if ((M(1) < (windowsize(3)-figattrib.sidePanelWidth+3)) && (M(1) > (windowsize(3)-figattrib.sidePanelWidth-3)))
    figattrib.panelDrag = 1;
end
%if the mouse click was within the graph window
if ((M(1)>figloc(1))&&(M(1)<(figloc(1)+figloc(3)))&&(M(2)>figloc(2))&&(M(2)<(figloc(2)+figloc(4)))&&(clustattrib.nodata == 0))
    clicktype = get(hObject,'SelectionType');
    if strcmp(clicktype,'extend')
        graphattrib.shiftclick = 1;
    else
        graphattrib.shiftclick = 0;
    end
    if strcmp(clicktype,'alt')
        graphattrib.rightclick = 1;
    else
        graphattrib.rightclick = 0;
    end
    
    graphattrib.squarehighlighton = 0;
    viewrange = get(graphattrib.graphwindow,'UserData');
    %action depends on which tool is selected
    switch figattrib.toolselect
        case 1 %polygon drawing tool
            if ((graphattrib.polydraw == 0)&&(~graphattrib.rightclick)) %this is the first point
                
                %we need to make sure that the current filtered cluster
                %settings won't cause a dependency loop.  If it does, then
                %the action is undone.
                tmpdepend = clustattrib.dependencies;
                tmpdepend(clustattrib.currclust, graphattrib.blackedOutClusters) = true;
                numclust = size(figattrib.mixcolor,1);
                %the new dependency matrix taken to the power of the number of clusters
                %will give non-zero values for all nodes that have infinate chidren or
                %parents
                generationmatrix = (tmpdepend^numclust);
                cyclicCalc1 = ones(1,numclust)*generationmatrix; %which have infinate childen?
                cyclicCalc2 = generationmatrix*ones(numclust,1); %which have infinate parents?
                %the clusters that have both infinate parents and children are inside a
                %dependency loop.
                nonzeros = find((cyclicCalc1 > 0)&(cyclicCalc2' > 0));
                
                if ~isempty(nonzeros)
                    errhandle = errordlg(['With the current cluster filters, this will create a depenency loop between clusters: ',num2str(nonzeros),'. Action aborted.']);                  
                    return
                end
                
                if (findex == -1)
                    findex = addfilterindex(filterline);
                end
                
                %if this box already exists, delete it (for all activiated
                %time filters)
                for i = 1:length(findex2)
                    try 
                        delete(graphattrib.polyg(a1,a2,findex2(i)).lines{clustattrib.currclust});
                    end
                    %turn highlighting off
                    try 
                        delete(graphattrib.polyg(a1,a2,findex2(i)).highlight{clustattrib.currclust});             
                    end
                    %remove the record of a box in these axes if one exists
                    try
                        clustattrib.clusters{clustattrib.currclust}.defineaxes = setdiff(clustattrib.clusters{clustattrib.currclust}.defineaxes,[a1 a2 findex2(i)],'rows');
                    catch
                        clustattrib.clusters{clustattrib.currclust}.defineaxes = [];
                    end
                    
                end
                clustattrib.clusters{clustattrib.currclust}.defineaxes = unique([clustattrib.clusters{clustattrib.currclust}.defineaxes;[a1 a2 findex]],'rows');
                
                %turn cluster on
                turnoncluster(clustattrib.currclust,handles);
                graphattrib.polydraw = 1;
                temp = get(graphattrib.graphwindow,'CurrentPoint');
               
                %'vertices' are the points on the current plot.
                %'relativevertices' are between 0 and 1, relative to the
                %total data range.
                %'lines' are the handles to the drawn lines.
                %'type' is either 1 or 2 for polygon or square.
                

                
                graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust} = repmat(temp(1,1:2),3,1);
                graphattrib.polyg(a1,a2,findex).relativevertices{clustattrib.currclust} = [(graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(:,1)+viewrange(3))/viewrange(1) (graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(:,2)+viewrange(5))/viewrange(2)];
                graphattrib.polyg(a1,a2,findex).lines{clustattrib.currclust} = line(graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(:,1),graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(:,2));
                set(graphattrib.polyg(a1,a2,findex).lines{clustattrib.currclust},'ButtonDownFcn','matclust(''polygon_ButtonDownFcn'',gcbo,guidata(gcbo))');
                set(graphattrib.polyg(a1,a2,findex).lines{clustattrib.currclust},'UserData',[a1 a2 findex clustattrib.currclust]);
                set(graphattrib.polyg(a1,a2,findex).lines{clustattrib.currclust},'Color',colorclip(figattrib.mixcolor(clustattrib.currclust,:)-[.2 .2 .2]));                
                set(graphattrib.polyg(a1,a2,findex).lines{clustattrib.currclust},'LineStyle','-','LineWidth',1);                

                %set(graphattrib.polyg(a1,a2,findex).lines{clustattrib.currclust},'UIContextMenu',handles.polycontext);
                graphattrib.polyg(a1,a2,findex).type{clustattrib.currclust} = 1;
                clustattrib.clusters{clustattrib.currclust}.box{a1,a2} = graphattrib.polyg(a1,a2,findex).relativevertices{clustattrib.currclust};    
            elseif (graphattrib.polydraw == 1) %we are in the middle of drawing a polygon
                if (stopwatch < .3) %double click ends the polygon drawing (here we define the sensitivity of the double-click in seconds)                                     
                 %if  strcmp(clicktype,'open')   %double click ends the drawing
                 
                    graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(end-1,1:2) = graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(end,1:2);
                    graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust} = graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(1:end-1,1:2);
                    graphattrib.polyg(a1,a2,findex).relativevertices{clustattrib.currclust} = [(graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(:,1)+viewrange(3))/viewrange(1) (graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(:,2)+viewrange(5))/viewrange(2)];                 
                    set(graphattrib.polyg(a1,a2,findex).lines{clustattrib.currclust},'UIContextMenu',handles.polycontext);
                    releasepolydraw(handles);                 
                else %not a double-click, so add new point to polygon
                    temp = get(graphattrib.graphwindow,'CurrentPoint');
                    graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(end-1,1:2) = temp(1,1:2);
                    graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(end,1:2) = temp(1,1:2);
                    graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(end+1,1:2) = graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(1,1:2);
                    graphattrib.polyg(a1,a2,findex).relativevertices{clustattrib.currclust} = [(graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(:,1)+viewrange(3))/viewrange(1) (graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(:,2)+viewrange(5))/viewrange(2)];                 
                    set(graphattrib.polyg(a1,a2,findex).lines{clustattrib.currclust},'XData',graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(:,1));
                    set(graphattrib.polyg(a1,a2,findex).lines{clustattrib.currclust},'YData',graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(:,2));                              
                    clustattrib.clusters{clustattrib.currclust}.box{a1,a2,findex} = graphattrib.polyg(a1,a2,findex).relativevertices{clustattrib.currclust};        
                    if strcmp(clicktype,'extend') %shift-click or middle click also end the polygon drawing, and the last point is included 
                        set(graphattrib.polyg(a1,a2,findex).lines{clustattrib.currclust},'UIContextMenu',handles.polycontext);
                        releasepolydraw(handles);
                    end
                end                               
            end
            
        case 2 %square drawing tool
            
            if (~graphattrib.rightclick)
                %we need to make sure that the current filtered cluster
                %settings won't cause a dependency loop.  If it does, then
                %the action is undone.
                tmpdepend = clustattrib.dependencies;
                tmpdepend(clustattrib.currclust, graphattrib.blackedOutClusters) = true;
                numclust = size(figattrib.mixcolor,1);
                %the new dependency matrix taken to the power of the number of clusters
                %will give non-zero values for all nodes that have infinate chidren or
                %parents
                generationmatrix = (tmpdepend^numclust);
                cyclicCalc1 = ones(1,numclust)*generationmatrix; %which have infinate childen?
                cyclicCalc2 = generationmatrix*ones(numclust,1); %which have infinate parents?
                %the clusters that have both infinate parents and children are inside a
                %dependency loop.
                nonzeros = find((cyclicCalc1 > 0)&(cyclicCalc2' > 0));
                
                if ~isempty(nonzeros)
                    errhandle = errordlg(['With the current cluster filters, this will create a depenency loop between clusters: ',num2str(nonzeros),'. Action aborted.']);
                    return
                end
                %the square is set up here, and made permanent once the mouse
                %button is released
                graphattrib.drawsquare = 1;
                temp = get(graphattrib.graphwindow,'CurrentPoint');
                if (findex == -1)
                    findex = addfilterindex(filterline);
                end
                for i = 1:length(findex2)
                    try
                        delete(graphattrib.polyg(a1,a2,findex2(i)).lines{clustattrib.currclust});
                    end
                    %turn highlighting off
                    try
                        delete(graphattrib.polyg(a1,a2,findex2(i)).highlight{clustattrib.currclust});
                    end
                    %remove the record of a box in these axes if one exists
                    try
                        clustattrib.clusters{clustattrib.currclust}.defineaxes = setdiff(clustattrib.clusters{clustattrib.currclust}.defineaxes,[a1 a2 findex2(i)],'rows');
                    catch
                        clustattrib.clusters{clustattrib.currclust}.defineaxes = [];
                    end
                    
                end
                
                %turn cluster on
                turnoncluster(clustattrib.currclust,handles);
                clustattrib.clusters{clustattrib.currclust}.defineaxes = unique([clustattrib.clusters{clustattrib.currclust}.defineaxes;[a1 a2 findex]],'rows');
                updateclustinfo(handles);
                graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust} = repmat(temp(1,1:2),5,1);
                graphattrib.polyg(a1,a2,findex).relativevertices{clustattrib.currclust} = [(graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(:,1)+viewrange(3))/viewrange(1) (graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(:,2)+viewrange(5))/viewrange(2)];
                graphattrib.polyg(a1,a2,findex).lines{clustattrib.currclust} = line(graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(1:4,1),graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(1:4,2));
                set(graphattrib.polyg(a1,a2,findex).lines{clustattrib.currclust},'ButtonDownFcn','matclust(''polygon_ButtonDownFcn'',gcbo,guidata(gcbo))');
                set(graphattrib.polyg(a1,a2,findex).lines{clustattrib.currclust},'Color',colorclip(figattrib.mixcolor(clustattrib.currclust,:)-[.2 .2 .2]));
                set(graphattrib.polyg(a1,a2,findex).lines{clustattrib.currclust},'UserData',[a1 a2 findex clustattrib.currclust]);
                set(graphattrib.polyg(a1,a2,findex).lines{clustattrib.currclust},'UIContextMenu',handles.polycontext);
                graphattrib.polyg(a1,a2,findex).type{clustattrib.currclust} = 1;
                clustattrib.clusters{clustattrib.currclust}.box{a1,a2,findex} = graphattrib.polyg(a1,a2,findex).relativevertices{clustattrib.currclust};
            end
        case 3 %arrow tool is used to select clustattrib.clusters
            %this section is called if the arrow is clicked in the graph. All highlighting is turned off
            %when this happens, except if the arrow hit a highlighted
            %square.
                     
            htest = get(gco,'UserData');
            try
                highlighthit = htest.type;
            catch
                highlighthit = 0;
            end          
            
            if (~highlighthit) %if a highlighted vertex was not hit, delete all highlighting                   
                
                clearhighlight;
             end
        case 4 %magnifyier tool
            %the magnifyer box is almost identical in code to the
            %cluter box above
            graphattrib.magdrawsquare = 1;
            temp = get(graphattrib.graphwindow,'CurrentPoint');    
            viewrange = get(graphattrib.graphwindow,'UserData');
            try 
                delete(graphattrib.magsquare.lines);
            end
                           
            graphattrib.magsquare.vertices = repmat(temp(1,1:2),5,1);
            graphattrib.magsquare.relativevertices = [((graphattrib.magsquare.vertices(:,1)+viewrange(3))/viewrange(1)) ((graphattrib.magsquare.vertices(:,2)+viewrange(5))/viewrange(2))];
            
            graphattrib.magsquare.lines = line(graphattrib.magsquare.vertices(1:4,1),graphattrib.magsquare.vertices(1:4,2));
            set(graphattrib.magsquare.lines,'Color',1-graphattrib.backgroundcolor);
            set(graphattrib.magsquare.lines,'LineStyle','--');
        case 5 %hand tool
            clearhighlight;
            graphattrib.plothighlightclear = 1;
            graphattrib.handdown = 1;
            temp = get(graphattrib.graphwindow,'CurrentPoint');    
            viewrange = get(graphattrib.graphwindow,'UserData');
            
            graphattrib.handinfo = [temp(1,1:2) viewrange viewrange(3:6)]; % orig cursor point -- orig viewrange -- new viewrange
            
        case 6 %magic wand
            if (~graphattrib.rightclick)
                %we need to make sure that the current filtered cluster
                %settings won't cause a dependency loop.  If it does, then
                %the action is undone.
                tmpdepend = clustattrib.dependencies;
                tmpdepend(clustattrib.currclust, graphattrib.blackedOutClusters) = true;
                numclust = size(figattrib.mixcolor,1);
                %the new dependency matrix taken to the power of the number of clusters
                %will give non-zero values for all nodes that have infinate chidren or
                %parents
                generationmatrix = (tmpdepend^numclust);
                cyclicCalc1 = ones(1,numclust)*generationmatrix; %which have infinate childen?
                cyclicCalc2 = generationmatrix*ones(numclust,1); %which have infinate parents?
                %the clusters that have both infinate parents and children are inside a
                %dependency loop.
                nonzeros = find((cyclicCalc1 > 0)&(cyclicCalc2' > 0));
                
                if ~isempty(nonzeros)
                    errhandle = errordlg(['With the current cluster filters, this will create a depenency loop between clusters: ',num2str(nonzeros),'. Action aborted.']);
                    return
                end
                
                if (findex == -1)
                    findex = addfilterindex(filterline);
                end
                
                %if this box already exists, delete it (for all activiated
                %time filters)
                for i = 1:length(findex2)
                    try
                        delete(graphattrib.polyg(a1,a2,findex2(i)).lines{clustattrib.currclust});
                    end
                    %turn highlighting off
                    try
                        delete(graphattrib.polyg(a1,a2,findex2(i)).highlight{clustattrib.currclust});
                    end
                    %remove the record of a box in these axes if one exists
                    try
                        clustattrib.clusters{clustattrib.currclust}.defineaxes = setdiff(clustattrib.clusters{clustattrib.currclust}.defineaxes,[a1 a2 findex2(i)],'rows');
                    catch
                        clustattrib.clusters{clustattrib.currclust}.defineaxes = [];
                    end
                    
                end
                clustattrib.clusters{clustattrib.currclust}.defineaxes = unique([clustattrib.clusters{clustattrib.currclust}.defineaxes;[a1 a2 findex]],'rows');
                
                
                set(figattrib.handles.statusText,'String','Creating polygon...');
                set(figattrib.handles.statusText,'UserData','Creating polygon...');
                drawnow;
                
                %turn cluster on
                turnoncluster(clustattrib.currclust,handles);
                
                temp = get(graphattrib.graphwindow,'CurrentPoint');
                relativePoint = [(temp(1,1)+viewrange(3))/viewrange(1) (temp(1,2)+viewrange(5))/viewrange(2)];
                
                
                filts = clustattrib.filterindex{findex};
                tmpfilter = fastorbit(clustdata.timefilters,filts(1)); %temporal filters
                tmpfilter2 = fastandbit(clustdata.otherfilters,filts(2:end)); %other filters
                
                if (~clustattrib.cluster0attrib.show)
                    visiblePoints = false(length(tmpfilter),1);
                    for i = 1:length(clustattrib.clusters)
                        if (~clustattrib.hiddenclusters(i,2))
                            visiblePoints(clustattrib.clusters{i}.index) = true;
                        end
                    end
                else
                    visiblePoints = true(length(tmpfilter),1);
                end
                
                %filteredpoints = find(tmpfilter & tmpfilter2 & graphattrib.nonhiddenpoints & visiblePoints);
                filteredpoints = find(clustdata.filteredpoints & graphattrib.nonhiddenpoints & visiblePoints);
                
                
                if (a1 > clustdata.origparams)
                    tmprotdata = clustdata.params(filteredpoints,clustdata.rotation(a1).params);
                    tmprotdata(:,4) = 1;
                    
                    
                    projection = (clustdata.rotation(a1).matrix*(tmprotdata'))';
                    plotspace = projection(:,[1 2]);
                    currentdatarange = clustdata.rotation(a1).datarange(1:2,[2 1]);
                    
                    clear tmprotdata;
                    clear projection;
                else
                    
                    plotspace = clustdata.params(filteredpoints,[a1 a2]);
                    currentdatarange = [clustdata.datarange(1,a2) clustdata.datarange(1,a1);clustdata.datarange(2,a2) clustdata.datarange(2,a1)];
                end
                
                
                xpoints = ((plotspace(:,1)-currentdatarange(1,2)))/(currentdatarange(2,2)-currentdatarange(1,2));
                ypoints = ((plotspace(:,2)-currentdatarange(1,1)))/(currentdatarange(2,1)-currentdatarange(1,1));
                               
                
                %this is where the polygon vertices are automatically computed
                
                
%                 try
%                     currPointsInCluster = clustattrib.clusters{clustattrib.currclust}.index;
%                 catch
%                     currPointsInCluster = [];
%                 end           
%                 definedInCurrentFilter = 0;                
%                 if sum(ismember(currPointsInCluster,filteredpoints))
%                     definedInCurrentFilter = 1;
%                 end
                
                dPoints = get(handles.wandSlider, 'Value');
                stepFrac = 5;
                if (clustattrib.cluster0attrib.show)
                    relativeVerts = findCloudBorders([xpoints ypoints],relativePoint, dPoints, stepFrac);
                else
                    %if cluster 0 is hidden, maybe run a different
                    %algorithm?
                    relativeVerts = findCloudBorders([xpoints ypoints],relativePoint, dPoints, stepFrac);
                end
                
                %'vertices' are the points on the current plot.
                %'relativevertices' are between 0 and 1, relative to the
                %total data range.
                %'lines' are the handles to the drawn lines.
                %'type' is either 1 or 2 for polygon or square.
                
                graphattrib.polyg(a1,a2,findex).relativevertices{clustattrib.currclust} = relativeVerts;
                graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust} = [(relativeVerts(:,1)*viewrange(1))-viewrange(3) (relativeVerts(:,2)*viewrange(2))-viewrange(5)];
                graphattrib.polyg(a1,a2,findex).lines{clustattrib.currclust} = line(graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(:,1),graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(:,2));
                set(graphattrib.polyg(a1,a2,findex).lines{clustattrib.currclust},'ButtonDownFcn','matclust(''polygon_ButtonDownFcn'',gcbo,guidata(gcbo))');
                set(graphattrib.polyg(a1,a2,findex).lines{clustattrib.currclust},'UserData',[a1 a2 findex clustattrib.currclust]);
                set(graphattrib.polyg(a1,a2,findex).lines{clustattrib.currclust},'Color',colorclip(figattrib.mixcolor(clustattrib.currclust,:)-[.2 .2 .2]));
                set(graphattrib.polyg(a1,a2,findex).lines{clustattrib.currclust},'LineStyle','-','LineWidth',1);
                set(graphattrib.polyg(a1,a2,findex).lines{clustattrib.currclust},'UIContextMenu',handles.polycontext);
                
                %set(graphattrib.polyg(a1,a2,findex).lines{clustattrib.currclust},'UIContextMenu',handles.polycontext);
                graphattrib.polyg(a1,a2,findex).type{clustattrib.currclust} = 1;
                clustattrib.clusters{clustattrib.currclust}.box{a1,a2} = graphattrib.polyg(a1,a2,findex).relativevertices{clustattrib.currclust};
                
                findInPoints(a1,a2,findex,clustattrib.currclust);
                if (length(find(clustdata.timefiltersOn)) > 1)
                    tfilts = find(clustdata.timefiltersOn);
                    for tnum = 2:length(tfilts)
                        copypoly(a1,a2,findex,clustattrib.currclust,tfilts(tnum)); %copy the polygon
                    end
                end
                updateclustinfo(handles);
                plotgraph(handles);
                highlightCurrentCluster(handles);
                addnewstate('create polygon',handles);
            end
            
    end
end
% -----------------------------------------------------------------
function figure1_WindowButtonMotionFcn(hObject, handles)
%called whenever the mouse moves within the figure window

global figattrib;
global graphattrib;
global clustattrib;
global clustdata;


a1 = get(handles.listbox1,'value');
a2 = get(handles.listbox2,'value');
[findex, filterline] = findfilterindex(handles);
findex = findex(1);
filterline = filterline(1,:);
[findex2, filterline2] = findfilterindex(handles);
figloc = get(graphattrib.graphwindow,'Position');
windowsize = get(handles.figure1,'Position');
M = get(hObject,'CurrentPoint');


if (figattrib.panelDrag) %we are dragging the side panel, so resize it as the cursor moves
    setPanelWidth(windowsize(3)-M(1),0);
end
    
if ((M(1)>figloc(1))&&(M(1)<(figloc(1)+figloc(3)))&&(M(2)>figloc(2))&&(M(2)<(figloc(2)+figloc(4)))&&(clustattrib.nodata == 0)) %is the motion within the plot window?
    cf = get(handles.figure1,'CurrentObject');
    cfprops = get(cf);
    try
       isuicontrol = strcmp(cfprops.Type,'uicontrol');
    catch
       isuicontrol = 0;
    end
    
    if (~figattrib.mainFigFocus && isuicontrol)
        %ReleaseFocus(handles.figure1); %annoying work-around to give main figure focus for keyboard shortcuts
        set(handles.figure1,'CurrentObject',handles.figure1);  
        figattrib.mainFigFocus = 1;
        
    end
    %Make the status bar display the x and y info of the cursor
    viewrange = get(graphattrib.graphwindow,'UserData');
    figpoint = get(graphattrib.graphwindow,'CurrentPoint');
    figpoint = figpoint(1,1:2);
    currviewbox = graphattrib.viewbox(a1,:,a2);
    s = get(graphattrib.graphwindow,'Position');
	resolutionx = s(3)*graphattrib.resolutionfactor(1);
	resolutiony = s(4)*graphattrib.resolutionfactor(2);
	currviewbox = graphattrib.viewbox(a1,:,a2);
	
	xrange = ((abs(diff(graphattrib.currentdatarange(:,2))))*[currviewbox(1);currviewbox(2)])+ (graphattrib.currentdatarange(1,2));
	yrange = ((abs(diff(graphattrib.currentdatarange(:,1))))*[currviewbox(3);currviewbox(4)])+ (graphattrib.currentdatarange(1,1));

    xpercentages = (figpoint(1)-1)/(resolutionx-1);
	ypercentages = (figpoint(2)-1)/(resolutiony-1);
	newx = round((((abs(diff(xrange)))*xpercentages)+xrange(1))*10)/10;
	if (a1 == 1)
        newx = timetrans(newx',clustdata.UnitsPerSec,1)';
        newx = newx{1};    
    else
        newx = num2str(newx);
    end
	newy = round((((abs(diff(yrange)))*ypercentages)+yrange(1))*10)/10;
	if (a2 == 1)
        newy = timetrans(newy',clustdata.UnitsPerSec,1)';
        newy = newy{1};    
    else
        newy = num2str(newy);
    end
    postext = [newx,' , ',newy];
    currstatus = get(handles.statusText,'UserData');
    %currstatus(12:12+length(postext)-1) = postext;
    currstatus = [currstatus,'       ',postext];
    
    set(handles.statusText,'String',currstatus);
   
    %Depending on which tool is currently selected, different things must
    %happen...
    switch figattrib.toolselect
      case 1 %polygon drawing tool
          set(hObject,'Pointer','crosshair');
          figpoint = get(graphattrib.graphwindow,'CurrentPoint');
          figpoint = figpoint(1,1:2);
          if (graphattrib.polydraw) %polygon currently being drawn
              graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(end-1,:) = figpoint;
              graphattrib.polyg(a1,a2,findex).relativevertices{clustattrib.currclust} = [(graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(:,1)+viewrange(3))/viewrange(1) (graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(:,2)+viewrange(5))/viewrange(2)];

              tmp = graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust};
              set(graphattrib.polyg(a1,a2,findex).lines{clustattrib.currclust},'XData',tmp(:,1));
              set(graphattrib.polyg(a1,a2,findex).lines{clustattrib.currclust},'YData',tmp(:,2));
              clustattrib.clusters{clustattrib.currclust}.box{a1,a2,findex} = graphattrib.polyg(a1,a2,findex).relativevertices{clustattrib.currclust};    
          end
      case 2 %square drawing tool
          set(hObject,'Pointer','crosshair');
          figpoint = get(graphattrib.graphwindow,'CurrentPoint');
          figpoint = figpoint(1,1:2);
          if (graphattrib.drawsquare) %square currently being drawn
              graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(2,1) = figpoint(1);
              graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(3,:) = figpoint;
              graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(4,2) = figpoint(2);
              graphattrib.polyg(a1,a2,findex).relativevertices{clustattrib.currclust} = [(graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(:,1)+viewrange(3))/viewrange(1) (graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(:,2)+viewrange(5))/viewrange(2)];
              set(graphattrib.polyg(a1,a2,findex).lines{clustattrib.currclust},'XData',graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(:,1));
              set(graphattrib.polyg(a1,a2,findex).lines{clustattrib.currclust},'YData',graphattrib.polyg(a1,a2,findex).vertices{clustattrib.currclust}(:,2));
              clustattrib.clusters{clustattrib.currclust}.box{a1,a2,findex} = graphattrib.polyg(a1,a2,findex).relativevertices{clustattrib.currclust};        
          end
      case 3 %arrow (selector) tool
          totalhighlightindex = [];
          set(hObject,'Pointer','arrow');
          figpoint = get(graphattrib.graphwindow,'CurrentPoint');
          figpoint = figpoint(1,1:2);
          if (graphattrib.squarehighlighton) %if the mouse is currently pressing on a polygons vertex (move vertex)
              totalhighlightindex = get(graphattrib.currentpolyhighlight,'UserData');
              if ~iscell(totalhighlightindex)
                  tmptotalhighlightindex{1} = totalhighlightindex;
                  totalhighlightindex = tmptotalhighlightindex;
              end
              for i = 1:length(totalhighlightindex)
                  highlightindex = totalhighlightindex{i};
                  findex = highlightindex(3); %make the filter index for the polygon instead of the current filter 
                  try
                    vertexindex = get(graphattrib.polyg(highlightindex(1),highlightindex(2),highlightindex(3)).highlight{highlightindex(4)},'UserData');
                  catch
                    continue;
                  end
                  vertexindex = vertexindex.index;
                  tmp = graphattrib.polyg(a1,a2,findex).vertices{highlightindex(4)};
                  
                  if ((vertexindex == size(tmp,1))|(vertexindex == 1)) %because the first and last points of the polygons are the same, they must be moved together
                      graphattrib.polyg(a1,a2,findex).vertices{highlightindex(4)}(1,1:2) = figpoint;
                      graphattrib.polyg(a1,a2,findex).vertices{highlightindex(4)}(size(tmp,1),1:2) = figpoint;
                  else %otherwise just move the single vertices
                      graphattrib.polyg(a1,a2,findex).vertices{highlightindex(4)}(vertexindex,1:2) = figpoint;
                  end
                  
                  graphattrib.polyg(a1,a2,findex).relativevertices{highlightindex(4)} = [(graphattrib.polyg(a1,a2,findex).vertices{highlightindex(4)}(:,1)+viewrange(3))/viewrange(1) (graphattrib.polyg(a1,a2,findex).vertices{highlightindex(4)}(:,2)+viewrange(5))/viewrange(2)];
                  graphattrib.pointmoved = highlightindex(4);
                  tmp = graphattrib.polyg(a1,a2,findex).vertices{highlightindex(4)};             
                  set(graphattrib.polyg(a1,a2,findex).lines{highlightindex(4)},'XData',tmp(:,1));
                  set(graphattrib.polyg(a1,a2,findex).lines{highlightindex(4)},'YData',tmp(:,2));
                  set(graphattrib.polyg(a1,a2,findex).highlight{highlightindex(4)},'Position',[tmp(vertexindex,1)-(3*graphattrib.resolutionfactor(1)) tmp(vertexindex,2)-(4*graphattrib.resolutionfactor(2)) 7*graphattrib.resolutionfactor(1) 7*graphattrib.resolutionfactor(2)]);
                  clustattrib.clusters{highlightindex(4)}.box{a1,a2,findex} = graphattrib.polyg(a1,a2,findex).relativevertices{highlightindex(4)}; 
              end
          end        
          if (graphattrib.polypress) %if the mouse is currently pressing on a polygon (not on a vertex), drag entire polygon
              totalhighlightindex = get(graphattrib.currentpolyhighlight,'UserData');
              if ~iscell(totalhighlightindex)
                  tmptotalhighlightindex{1} = totalhighlightindex;
                  totalhighlightindex = tmptotalhighlightindex;
              end
              for i = 1:length(totalhighlightindex)
                  highlightindex = totalhighlightindex{i};
                  
                  findex = highlightindex(3);
                  tmp = graphattrib.polyg(a1,a2,findex).vertices{highlightindex(4)};
                  graphattrib.polyg(a1,a2,findex).vertices{highlightindex(4)} = [figpoint(1)-graphattrib.relativepolypress{i}(:,1) figpoint(2)-graphattrib.relativepolypress{i}(:,2)];
                  graphattrib.polyg(a1,a2,findex).relativevertices{highlightindex(4)} = [(graphattrib.polyg(a1,a2,findex).vertices{highlightindex(4)}(:,1)+viewrange(3))/viewrange(1) (graphattrib.polyg(a1,a2,findex).vertices{highlightindex(4)}(:,2)+viewrange(5))/viewrange(2)];
                  
                  graphattrib.pointmoved = highlightindex(4);
                  tmp = graphattrib.polyg(a1,a2,findex).vertices{highlightindex(4)};
                  
                  set(graphattrib.polyg(a1,a2,findex).lines{highlightindex(4)},'XData',tmp(:,1));
                  set(graphattrib.polyg(a1,a2,findex).lines{highlightindex(4)},'YData',tmp(:,2));
                  clustattrib.clusters{highlightindex(4)}.box{a1,a2,findex} = graphattrib.polyg(a1,a2,findex).relativevertices{highlightindex(4)}; 
              end
          end   

      case 4 %magnifying tool
          %set(hObject,'Pointer','circle');
          load magpic;
          set(hObject,'Pointer','custom','PointerShapeCData',magpic,...
            'PointerShapeHotSpot',[6 5]);
          figpoint = get(graphattrib.graphwindow,'CurrentPoint');
          figpoint = figpoint(1,1:2);
          if (graphattrib.magdrawsquare) %if a magifying square is currently being drawn 
              graphattrib.magsquare.vertices(2,1) = figpoint(1);
              graphattrib.magsquare.vertices(3,:) = figpoint;
              graphattrib.magsquare.vertices(4,2) = figpoint(2);
              graphattrib.magsquare.relativevertices = [((graphattrib.magsquare.vertices(:,1)+viewrange(3))/viewrange(1)) ((graphattrib.magsquare.vertices(:,2)+viewrange(5))/viewrange(2))];

              set(graphattrib.magsquare.lines,'XData',graphattrib.magsquare.vertices(:,1));
              set(graphattrib.magsquare.lines,'YData',graphattrib.magsquare.vertices(:,2));
          end
      case 5
          %set(hObject,'Pointer','circle');
          load handpointer;
          set(hObject,'Pointer','custom','PointerShapeCData',handpointer,...
            'PointerShapeHotSpot',[8 8]);
         figpoint = get(graphattrib.graphwindow,'CurrentPoint');
         figpoint = figpoint(1,1:2);
         if (graphattrib.handdown)
            s = get(graphattrib.graphwindow,'Position');
            resolutionx = s(3)*graphattrib.resolutionfactor(1);
            resolutiony = s(4)*graphattrib.resolutionfactor(2);

            spplotspace = figattrib.spplotspace; 
            
            change = figpoint - graphattrib.handinfo(1:2);
            change2 = [change(1) change(1) change(2) change(2)];
            
            viewrange = graphattrib.handinfo(5:8)-change2;
           %viewrange(2) = min([viewrange(2) size(spplotspace,2)]);
           %viewrange(4) = min([viewrange(4) size(spplotspace,1)]);
           %size(spplotspace)
               
           %size(spplotspace)
           viewrange = round(viewrange);
           graphattrib.handinfo(3:4);
            if ((viewrange(1) >= 1) & (viewrange(2) <=  size(spplotspace,2)))
                %viewrange = [max([1 viewrange(1)]) min([size(spplotspace,2) viewrange(2)]) max([1 viewrange(3)]) min([size(spplotspace,1) viewrange(4)]) ];
                    set(graphattrib.Himage,'CData',spplotspace(graphattrib.handinfo(11):graphattrib.handinfo(12),viewrange(1):viewrange(2)));
                    graphattrib.handinfo(9:10) = viewrange(1:2);
            end

            if ((viewrange(3) >= 1) & (viewrange(4) <=  size(spplotspace,1)))
                %viewrange = [max([1 viewrange(1)]) min([size(spplotspace,2) viewrange(2)]) max([1 viewrange(3)]) min([size(spplotspace,1) viewrange(4)]) ];
                    set(graphattrib.Himage,'CData',spplotspace(viewrange(3):viewrange(4),graphattrib.handinfo(9):graphattrib.handinfo(10)));
                    graphattrib.handinfo(11:12) = viewrange(3:4);              
            end
             %axis([1 resolutionx 1 resolutiony])
            currrange = [graphattrib.handinfo(3:4) graphattrib.handinfo(9:12)]; %current viewrange
            for tnum = 1:length(findex2)
                try
                    for i = 1:length(graphattrib.polyg(a1,a2,findex2(tnum)).vertices)
                        try
                            graphattrib.polyg(a1,a2,findex2(tnum)).vertices{i} = [(graphattrib.polyg(a1,a2,findex2(tnum)).relativevertices{i}(:,1)*currrange(1))-currrange(3) (graphattrib.polyg(a1,a2,findex2(tnum)).relativevertices{i}(:,2)*currrange(2))-currrange(5)];
                            set(graphattrib.polyg(a1,a2,findex2(tnum)).lines{i},'XData',graphattrib.polyg(a1,a2,findex2(tnum)).vertices{i}(:,1));
                            set(graphattrib.polyg(a1,a2,findex2(tnum)).lines{i},'YData',graphattrib.polyg(a1,a2,findex2(tnum)).vertices{i}(:,2));
                            set(graphattrib.polyg(a1,a2,findex2(tnum)).lines{i},'Color',figattrib.mixcolor(i,:));    
                        end
                    end
                end
            end
	
            
%              if ((viewrange(1) >= 1) & (viewrange(2) <  size(spplotspace,2)) & (viewrange(3) >= 1) & (viewrange(4) <  size(spplotspace,2)))
%                 %viewrange = [max([1 viewrange(1)]) min([size(spplotspace,2) viewrange(2)]) max([1 viewrange(3)]) min([size(spplotspace,1) viewrange(4)]) ];
%                 
%                 set(graphattrib.Himage,'CData',spplotspace(viewrange(3):viewrange(4),viewrange(1):viewrange(2)));
%                 graphattrib.handinfo(9:12) = viewrange;
%                 %axis([0 resolutionx 0 resolutiony])
%             end
%             
            
            
            %axis(graphattrib.handinfo(3:6)-change2);
            %graphattrib.viewbox(a1,1:4,a2)     
         end
    case 6 %wand tool
          
          load wandpic;
          set(hObject,'Pointer','custom','PointerShapeCData',wandpic,...
            'PointerShapeHotSpot',[6 6]);
          figpoint = get(graphattrib.graphwindow,'CurrentPoint');
          figpoint = figpoint(1,1:2);
           
    end
else
    if ((M(1) < (windowsize(3)-figattrib.sidePanelWidth+3)) && (M(1) > (windowsize(3)-figattrib.sidePanelWidth-3)))
        set(hObject,'Pointer','left');
    else
        set(hObject,'Pointer','arrow');
    end
    set(handles.statusText,'String',get(handles.statusText,'UserData'));
    figattrib.mainFigFocus = 0;
        
end
%------------------------------------------------------------------ 
function figure1_WindowButtonUpFcn(hObject, handles)
%called whenever the left mouse button is released inside the figure window

global graphattrib;
global clustattrib;
global clustdata;
global figattrib;


currx = get(handles.listbox1,'Value');
curry = get(handles.listbox2,'Value');
[findex, filterline] = findfilterindex(handles);
findex = findex(1);
filterline = filterline(1,:);
[findex2, filterline2] = findfilterindex(handles);
graphattrib.squarehighlighton = 0;
graphattrib.polypress = 0;
windowsize = get(handles.figure1,'Position');
M = get(hObject,'CurrentPoint');
if (figattrib.panelDrag) %we are finished dragging the side panel, save in defaults file
    figattrib.panelDrag = 0;
    setPanelWidth(windowsize(3)-M(1),1);   
end


if (graphattrib.drawsquare) %square was being drawn
    graphattrib.drawsquare = 0;
    findInPoints(currx,curry,findex,clustattrib.currclust); %redefine the cluster
    %if more than one time filter is selected, copy the
    %polygon
    if (length(find(clustdata.timefiltersOn)) > 1)
        tfilts = find(clustdata.timefiltersOn);
        for tnum = 2:length(tfilts)                        
            copypoly(currx,curry,findex,clustattrib.currclust,tfilts(tnum)); %copy the polygon
        end
    end
    updateclustinfo(handles);
    plotgraph(handles) %redraw the axes
    highlightCurrentCluster(handles);
    addnewstate('draw square',handles); 
end

if (graphattrib.magdrawsquare) %a magnifying square was being drawn
    graphattrib.magdrawsquare = 0;
    try
        delete(graphattrib.magsquare.lines);
    end
    tempbox = graphattrib.magsquare.relativevertices;
    graphattrib.oldviewbox = graphattrib.viewbox;
    graphattrib.viewbox(currx,1:4,curry) = [min(tempbox(:,1)) max(tempbox(:,1)) min(tempbox(:,2)) max(tempbox(:,2))];  %change the current view window to match the mag box    
    plotgraph(handles) %redraw axes
    addnewstate('zoom',handles); 
end   

if (graphattrib.pointmoved)  %if either a polygon vertex or an entire polygon has been moved during the mouse click
    
    %we need to make sure that the current filtered cluster
    %settings won't cause a dependency loop.  If it does, then
    %the action is undone.
    tmpdepend = clustattrib.dependencies;
    tmpdepend(graphattrib.pointmoved, graphattrib.blackedOutClusters) = true;
    numclust = size(figattrib.mixcolor,1);
    %the new dependency matrix taken to the power of the number of clusters
    %will give non-zero values for all nodes that have infinate chidren or
    %parents
    generationmatrix = (tmpdepend^numclust);
    cyclicCalc1 = ones(1,numclust)*generationmatrix; %which have infinate childen?
    cyclicCalc2 = generationmatrix*ones(numclust,1); %which have infinate parents?
    %the clusters that have both infinate parents and children are inside a
    %dependency loop.
    nonzeros = find((cyclicCalc1 > 0)&(cyclicCalc2' > 0));
    
    if ~isempty(nonzeros)
        errhandle = errordlg(['With the current cluster filters, this will create a depenency loop between clusters: ',num2str(nonzeros),'. Action aborted.']);
        waitfor(errhandle);
        clustattrib.currstate = clustattrib.currstate+1;
        undoMenu_Callback(0,handles);
        return
    end
    totalhighlightindex = get(graphattrib.currentpolyhighlight,'UserData');
    if ~iscell(totalhighlightindex)
        tmptotalhighlightindex{1} = totalhighlightindex;
        totalhighlightindex = tmptotalhighlightindex;
    end
    for i = 1:length(totalhighlightindex)
        highlightindex = totalhighlightindex{i};
        findex = highlightindex(3);
        
        %Because we are editing, and not re-creating a polygon,
        %we want to remember which clusters were filtered out at the time 
        %of drawing and keep these points out of the cluster.
        
        clustexcludes = find(clustattrib.dependencies(graphattrib.pointmoved,:));
        set(graphattrib.polyg(currx,curry,findex).lines{graphattrib.pointmoved},'LineStyle','-','LineWidth',1); 
        graphattrib.polyg(currx,curry,findex).type{graphattrib.pointmoved} = 1;
        findInPoints(currx,curry,findex,graphattrib.pointmoved,clustexcludes); %redefine the altered cluster
    end
    graphattrib.plothighlightclear = 0;
    plotgraph(handles) %redraw the axes
    graphattrib.pointmoved = 0;
    addnewstate('edit polygon',handles); 
end


if (graphattrib.handdown)
    graphattrib.handdown = 0;
    oldviewrange = get(graphattrib.graphwindow,'UserData');
    
    viewrange = [oldviewrange(1:2) graphattrib.handinfo(9:12)];
    relativerange = [viewrange(3)/viewrange(1) viewrange(4)/viewrange(1) viewrange(5)/viewrange(2) viewrange(6)/viewrange(2)];
    set(graphattrib.graphwindow,'UserData', viewrange);
    graphattrib.oldviewbox = graphattrib.viewbox;
    graphattrib.viewbox(currx,1:4,curry) = relativerange;  %change the current view window to match the mag box    
    currviewbox = relativerange;
    graphattrib.handinfo = [];
    s = get(graphattrib.graphwindow,'Position');
    resolutionx = s(3)*graphattrib.resolutionfactor(1);
    resolutiony = s(4)*graphattrib.resolutionfactor(2);
    XTick = get(graphattrib.graphwindow,'XTick');
	YTick = get(graphattrib.graphwindow,'YTick');
	newXTick = [1 resolutionx/4  resolutionx/2 3*resolutionx/4 resolutionx];
	newYTick = [1 resolutiony/4  resolutiony/2 3*resolutiony/4 resolutiony];
	set(graphattrib.graphwindow,'XTick',newXTick);
	set(graphattrib.graphwindow,'YTick',newYTick);
	XTick = newXTick;
	YTick = newYTick;
	
    xrange = ((abs(diff(graphattrib.currentdatarange(:,2))))*[currviewbox(1);currviewbox(2)])+ (graphattrib.currentdatarange(1,2));
	yrange = ((abs(diff(graphattrib.currentdatarange(:,1))))*[currviewbox(3);currviewbox(4)])+ (graphattrib.currentdatarange(1,1));
	xpercentages = (XTick-min(XTick))/(max(XTick)-min(XTick));
	ypercentages = (YTick-min(YTick))/(max(YTick)-min(YTick));
	newxlabels = round((((abs(diff(xrange)))*xpercentages)+xrange(1))*10)/10;
	if (currx == 1)
        newxlabels = timetrans(newxlabels',clustdata.UnitsPerSec,1)';
	end
	newylabels = round((((abs(diff(yrange)))*ypercentages)+yrange(1))*10)/10;
	if (curry == 1)
        newylabels = timetrans(newylabels',clustdata.UnitsPerSec,1)';
	end
	
	set(graphattrib.graphwindow,'XTickLabel',newxlabels);
	set(graphattrib.graphwindow,'YTickLabel',newylabels);
	set(graphattrib.graphwindow,'FontSize',7);
    
    %plotgraph(handles) %redraw the axes
    addnewstate('travel');
end 
%----------------------------------------------------------------
function polygon_ButtonDownFcn(hObject,handles)
%called whenever any of the polygons are clicked

global graphattrib;
global figattrib;
global clustattrib;

polyindex = get(hObject,'UserData');
currentpoly = graphattrib.polyg(polyindex(1),polyindex(2),polyindex(3)).vertices{polyindex(4)};
filterindex = polyindex(3);
objecthandle = graphattrib.polyg(polyindex(1),polyindex(2),polyindex(3)).lines{polyindex(4)};


if ~(graphattrib.shiftclick)
    [findex, filterline] = findfilterindex(handles);
    objecthandle = [];
    tempfilterindex = [];
    for i = 1:length(findex)
        try
            
            objecthandle = [objecthandle graphattrib.polyg(polyindex(1),polyindex(2),findex(i)).lines{polyindex(4)}];
            if ~isempty(graphattrib.polyg(polyindex(1),polyindex(2),findex(i)).lines{polyindex(4)})
                tempfilterindex = [tempfilterindex findex(i)];
            end
            
        end
    end
    filterindex = tempfilterindex;
end
%polyindex
%filterindex

figpoint = get(graphattrib.graphwindow,'CurrentPoint');
figpoint = figpoint(1,1:2);

if ((graphattrib.polydraw == 0)&&((figattrib.toolselect == 3)||(graphattrib.rightclick))) %arrow tool is selected
    graphattrib.relativepolypress = []; %location of each vertex relative to mouse (used while dragging the polygon)

    for i = 1:length(filterindex)
        currentpoly = graphattrib.polyg(polyindex(1),polyindex(2),filterindex(i)).vertices{polyindex(4)};
        try
            [dist(i),ind(i)] = min(sqrt((currentpoly(:,1)-figpoint(1)).^2 + (currentpoly(:,2)-figpoint(2)).^2)); %find distance from mouse to closest vertex
            graphattrib.relativepolypress{i} = [figpoint(1)-currentpoly(:,1) figpoint(2)-currentpoly(:,2)]; %location of each vertex relative to mouse (used while dragging the polygon)    
        catch
            dist(i) = 10000;
            ind(i) = 0;
        end
    end
    hold on
    if (polyindex(4) ~= clustattrib.currclust)
        figattrib.cashedPolygon = [];
        set(handles.pastePolygonFilterMainmenu,'Enable','off');
    end
    clustpickbutton_Callback(figattrib.clustcontrol(polyindex(4),2),handles)
    clearhighlight;
    
    
    
   
    
    if ~isempty(filterindex)
        graphattrib.polypress = 1; %polygon is being pressed
        vertexhits = [];
        for i = 1:length(filterindex)
            currentpoly = graphattrib.polyg(polyindex(1),polyindex(2),filterindex(i)).vertices{polyindex(4)};
            if (dist(i) < 10) %user clicked on a vertex, so highlight the vertex too
                vertexhits = [vertexhits i];
                hdata.type = 1;
                hdata.index = ind(i);
                hdata.polyindex = [polyindex(1) polyindex(2) filterindex(i) polyindex(4)];
                %the vertex highlighting is a filled square
                graphattrib.polyg(polyindex(1),polyindex(2),filterindex(i)).highlight{polyindex(4)} =  ...
                    rectangle('Position',[currentpoly(ind(i),1)-(3*graphattrib.resolutionfactor(1)) currentpoly(ind(i),2)-(4*graphattrib.resolutionfactor(2)) 7*(graphattrib.resolutionfactor(1)) 7*(graphattrib.resolutionfactor(2))],'FaceColor', figattrib.mixcolor(polyindex(4),:), ...
                    'EdgeColor',figattrib.mixcolor(polyindex(4),:),'ButtonDownFcn','matclust(''highlight_ButtonDownFcn'',gcbo,guidata(gcbo))', ...
                    'UserData',hdata);
                graphattrib.squarehighlighton = 1; %vertex is being pressed
                graphattrib.polypress = 0;  %non-vertex polygon is now NOT being pressed
            end
        end
        if (graphattrib.squarehighlighton)
            objecthandle = objecthandle(vertexhits);
        end
    end
    hold off
     
    set(objecthandle,'Marker','s'); %highlight the clicked polygon
    graphattrib.currentpolyhighlight = objecthandle; 
    
    if (length(objecthandle) == 1)
        %set(handles.mainPolygonMenu,'Enable','on');
        set(handles.polydeletemenu,'Enable','on');
        set(handles.addPolygonFilterMainmenu,'Enable','on');
        set(handles.copyPolygonFilterMainmenu,'Enable','on');

        Child = get(handles.polycontext,'Children');
        for ch = 1:length(Child)
            set(Child(ch),'Enable','on');
        end
    else
        %set(handles.mainPolygonMenu,'Enable','off');
        set(handles.polydeletemenu,'Enable','on');
        set(handles.addPolygonFilterMainmenu,'Enable','off');
        set(handles.copyPolygonFilterMainmenu,'Enable','off');
        Child = get(handles.polycontext,'Children');
        for ch = 1:length(Child)
            set(Child(ch),'Enable','off');
        end
        set(Child(2),'Enable','on');
    end
    
  

    %set(handles.mainClustMenu,'Enable','on');
    %for menuNum = 1:length(handles.mainClusterMenu)
    %    set(handles.mainClusterMenu(menuNum),'UserData',polyindex(4))
    %end
    
end
%----------------------------------------------------------
function figure1_WindowScrollWheelFcn(hObject,event)

global figattrib;
global graphattrib;
global clustattrib;
global clustdata;

handles = guidata(hObject);
figloc = get(graphattrib.graphwindow,'Position');
windowsize = get(handles.figure1,'Position');
M = get(hObject,'CurrentPoint');
    
if ((M(1)>figloc(1))&&(M(1)<(figloc(1)+figloc(3)))&&(M(2)>figloc(2))&&(M(2)<(figloc(2)+figloc(4)))&&(clustattrib.nodata == 0)) %is the cursor within the plot window?

    if (figattrib.shiftpress)
       
        sliderVal = get(handles.wandSlider,'value');
        sliderStep = get(handles.wandSlider,'SliderStep');
        sliderMax = get(handles.wandSlider,'Max');
        sliderMin = get(handles.wandSlider,'Min');
        
        sliderVal = sliderVal-event.VerticalScrollCount;
        if (sliderVal > sliderMax)
            sliderVal = sliderMax;
        end
        if (sliderVal < sliderMin)
            sliderVal = sliderMin;
        end
        set(handles.wandSlider,'value',sliderVal);
            
    else
        
        if (event.VerticalScrollCount < 1)
            zoomIn(M);
        else
            
            zoomOut(M);
        end
    end
end
%------------------------------------------------------
function cashPolygon(hObject, handles)

global graphattrib;
global figattrib;


highlightindex = get(graphattrib.currentpolyhighlight,'UserData');

%if (~isempty(highlightindex) && ~figattrib.tFiltAllowMultiple)
if (~isempty(highlightindex))
    
    figattrib.cashedPolygon = highlightindex;
    set(handles.pastePolygonFilterMainmenu,'Enable','on');
    
    
end
%------------------------------------------------------
function pasteCashedPolygon(hObject, handles)

global graphattrib;
global figattrib;
global clustdata;
global clustattrib;

if ((~isempty(figattrib.cashedPolygon)))
    
    [findex, filterline] = findfilterindex(handles);  
    tfilts = find(clustdata.timefiltersOn);
    for tnum = 1:length(tfilts)
                   
        if ~((findex(tnum) == figattrib.cashedPolygon(3)) && (clustattrib.currclust == figattrib.cashedPolygon(4)))
            try
                deletepoly(figattrib.cashedPolygon(1),figattrib.cashedPolygon(2),findex(tnum),clustattrib.currclust);
            end
            linepick = copypolyUniversal(figattrib.cashedPolygon(1),figattrib.cashedPolygon(2),figattrib.cashedPolygon(3),figattrib.cashedPolygon(4),clustattrib.currclust,filterline(tnum,:));
            %copypoly(figattrib.cashedPolygon(1),figattrib.cashedPolygon(2),figattrib.cashedPolygon(3),figattrib.cashedPolygon(4),tfilts(tnum));
        end
    end
    plotgraph(handles);
    highlightCurrentCluster(handles);
    FilterPoints;
    updateclustinfo(handles);
    addnewstate('copy polygon',handles);
    
end
%-------------------------------------------------------
function figure1_KeyReleaseFcn(hObject,event)

global figattrib;
global graphattrib;
global clustattrib;
handles = guidata(hObject);
if (figattrib.shiftpress)
    
    if ((~isempty(figattrib.cashedPolygon))&&(~figattrib.tFiltAllowMultiple)&&(strcmp(event.Key,'shift')))
        
        [findex, filterline] = findfilterindex(handles);
        if (findex ~= figattrib.cashedPolygon(3))
            try
                deletepoly(figattrib.cashedPolygon(1),figattrib.cashedPolygon(2),findex,figattrib.cashedPolygon(4));
            end
            linepick = copypolyUniversal(figattrib.cashedPolygon(1),figattrib.cashedPolygon(2),figattrib.cashedPolygon(3),figattrib.cashedPolygon(4),clustattrib.currclust,filterline);
            
            %copypoly(figattrib.cashedPolygon(1),figattrib.cashedPolygon(2),figattrib.cashedPolygon(3),figattrib.cashedPolygon(4),filterline(1))
            plotgraph(handles);
            objecthandle = graphattrib.polyg(figattrib.cashedPolygon(1),figattrib.cashedPolygon(2),linepick).lines{figattrib.cashedPolygon(4)};
            set(objecthandle,'Marker','s'); %highlight the clicked polygon
            graphattrib.currentpolyhighlight = objecthandle; 
            %set(handles.mainPolygonMenu,'Enable','on');
            set(handles.polydeletemenu,'Enable','on');
            set(handles.addPolygonFilterMainmenu,'Enable','on');
            set(handles.copyPolygonFilterMainmenu,'Enable','on');
            Child = get(handles.polycontext,'Children');
            for ch = 1:length(Child)
                set(Child(ch),'Enable','on');
            end
            FilterPoints;
            updateclustinfo(handles);
            addnewstate('copy polygon',handles); 
        end
        figattrib.cashedPolygon = [];
        set(handles.pastePolygonFilterMainmenu,'Enable','off');
        
    end
    if (strcmp(event.Key,'shift'))
        figattrib.shiftpress = 0;
    end
    
end



%-------------------------------------------------------
function figure1_KeyPressFcn(hObject,event)
%called when user hits a keyboard key


global clustdata;
global graphattrib;
global clustattrib;
global figattrib;

keypress = get(hObject,'CurrentCharacter');
numkeypress = double(keypress); %turn the character into an ASCII number
handles = guidata(hObject);
highlightindex = [];

if (length(event.Modifier) > 0)
    if strcmp(event.Modifier{1},'shift')
        figattrib.shiftpress = 1;
        try
            highlightindex = get(graphattrib.currentpolyhighlight,'UserData');
        end
        if (~isempty(highlightindex) && ~figattrib.tFiltAllowMultiple)

            figattrib.cashedPolygon = highlightindex;

        end
    end
end

controlPress = 0;
if (length(event.Modifier) > 0) && strcmp(event.Modifier{1},'control')
    controlPress = 1;
end

try
	switch numkeypress
        case {127, 8} %delete or backspace key
            if (~isempty(graphattrib.currentpolyhighlight)) %if a polygon is currently highlighted, delete it
                polyhandle = graphattrib.currentpolyhighlight;
                totalpolyindex = get(polyhandle,'UserData');
                if ~iscell(totalpolyindex)
                    tmppolyindex{1} = totalpolyindex;
                    totalpolyindex = tmppolyindex;
                end
                for i = 1:length(totalpolyindex)
                    polyindex = totalpolyindex{i};
                    deletepoly(polyindex(1),polyindex(2),polyindex(3),polyindex(4));
                end
                updateclustinfo(handles);
                plotgraph(handles);      
                addnewstate('delete polygon',handles); 
            end
        case 27 %escape key
            if (graphattrib.polydraw == 1)
                a1 = get(handles.listbox1,'value');
                a2 = get(handles.listbox2,'value');
                [findex, filterline] = findfilterindex(handles);
                try 
                    delete(graphattrib.polyg(a1,a2,findex).lines{clustattrib.currclust});
                end
                graphattrib.polydraw = 0;
                clustattrib.currstate = clustattrib.currstate+1;
                undoMenu_Callback(0,handles);

            end
        case 13
            if (graphattrib.polydraw == 1)
               releasepolydraw(handles);
            end
            
        %this case added by dylan 30h jan 2012    
        case {28,29}  %left/right key presses
            
            if (controlPress)
                %flip though the x,y combinations
                v1  = get(handles.listbox1,'value');
                v2  = get(handles.listbox2,'value');
                ex = numel(get(handles.listbox1,'string'));      %number of elements
                %ex = sum(cellfun(@(x)isempty(x),regexp(get(handles.listbox1,'string'),'ROTATION')))    %number of elements not containing the
                
                %define the possible combinations here
                %use a heuristic to work out if the user has put time as the first feature,
                %in which case it should be exluded from the list
                if any(diff(clustdata.params(1:(min([100 size(clustdata.params,1)])),1)) < 0)
                    cmbs = combnk([1:ex],2);
                else
                    cmbs = combnk([2:ex],2);
                end
                
                %exclude any combinations whos label includes the string
                %ROTATION and add them as separate pairs
                rot = find(cellfun(@(x)~isempty(x),regexp(get(handles.listbox1,'string'),'ROTATION')));
                exclude = sum(ismember(cmbs,rot),2) >0;
                cmbs(exclude,:) = [];
                cmbs = vertcat(cmbs,[rot rot]);
                
                
                %enumerate all possible permutation
                prms = [ceil((1:(ex*ex))./ex)' repmat((1:ex)',ex,1)];
                
                %find which one matches our current settings
                ind = find( prms(:,1) == v1 & prms(:,2) == v2);
                
                if numkeypress == 28;adj = -1;end   %left key
                if numkeypress == 29;adj = +1;end   %right key
                
                while true
                    %deal with overflow
                    ind = mod((ind-1)+adj,ex^2)+1;
                    
                    %try and find a match with the current prms and any of the
                    %cmbs (take the abs of the differences)
                    hit = sum(abs(repmat(prms(ind,:),size(cmbs,1),1) - cmbs),2)==0;
                    
                    %update the values and break
                    if any(hit)
                        v1 = prms(ind,1);
                        v2 = prms(ind,2);
                        break;
                    end
                end
                
                set(handles.listbox1,'value',v1);
                set(handles.listbox2,'value',v2);
                
                listbox1_Callback(handles.listbox1, [], handles)
                
                %if the last option we were in was a rotation we need to call
                %these again for the other list box, some pairing switch is
                %potentially not reset on these calls?
                if get(handles.listbox1,'value')~=v1 || get(handles.listbox2,'value')~=v2
                    set(handles.listbox1,'value',v1);
                    set(handles.listbox2,'value',v2);
                    listbox1_Callback(handles.listbox2, [], handles)
                end
            else
                
                %flip through the tools
                if numkeypress == 28;adj = -1;end   %left key
                if numkeypress == 29;adj = +1;end   %right key
                currentTool = figattrib.toolselect;
                newTool = mod(currentTool+adj-1,6)+1;
                figattrib.toolselect = newTool;
                switch newTool
                    case 1                     
                        polygonbutton_Callback();
                    case 2
                        squarebutton_Callback()
                    case 3
                        arrowbutton_Callback();
                    case 4
                        magbutton_Callback();
                    case 5
                        handbutton_Callback();
                    case 6
                        wandbutton_Callback();
                end
                
            end
            
            %this case added by dylan 30h jan 2012
        case {30,31}  %up down key presses 
            %switchesbetween the timefilters
            
            
            val = get(handles.timefilterList,'value');
            if isempty(val)
                val = 1;
            end
            
            nVal = size(clustdata.timefilterranges,1);
            if numkeypress == 30;adj = -1;end   %up
            if numkeypress == 31;adj = +1;end   %down
            
            val = val + adj;
            if (val > nVal); val = nVal; end
            if (val < 1);val = 1; end
            %val = mod((val-1+adj),nVal)+1
            set(handles.timefilterList,'value',val);            
            matclust('timefilterList_Callback',handles.timefilterList,handles)
            
            
          
            
        %this case added by dylan 4th FEB 2012    
        case {48,49,50,51,52,53,54,55,56,57}   %digits zero through 9
            if isempty(event.Modifier)
                %get the uicontrol handle for the cluster that refers to the
                %button hte user just pressed
                
                clusNum = numkeypress - 48;
                
                %now check that there was not a key pressed in the last few
                %half a second
                
                %             if isfield(figattrib,'lastKeyClustTog') && ...            %the field exists
                %                ~isempty(figattrib.lastKeyClustTog) && ...           %it is not empty
                %                toc(figattrib.lastKeyClustTog{1}) < 0.5              %it is hot on the tail of another key press
                %
                %                 %there was a key press in the last n seconds
                %                 %therefore the user wanted a two digit cluster
                %
                %                 %save the first key press
                %                 firstClust = figattrib.lastKeyClustTog{2};
                %
                %                 %generate the new two digit cluster number
                %                 clusNum = figattrib.lastKeyClustTog{2}*10 + clusNum;
                %
                %                 %clear the data since we are doiing a two digit cluster now
                %                 figattrib.lastKeyClustTog = [];
                %
                %                %recursively call the function so that we can toggle the
                %                %original cluster
                %                set(hObject,'CurrentCharacter',num2str( firstClust ));
                %                matclust('figure1_KeyPressFcn',hObject,handles);
                %
                %                %clear the data again since the last call would have filled
                %                %it up
                %                figattrib.lastKeyClustTog = [];
                %
                %
                %
                %             else
                %                 %store the data in case there is another character
                %                 figattrib.lastKeyClustTog{1} = tic;     %the ic here has its own scope and should not affet other tic timers
                %                 figattrib.lastKeyClustTog{2} = clusNum;
                %
                %             end
                
                if clusNum == 0
                    %toggle the state of the
                    newVal = (get(handles.cluster0button,'value')==0)+0;
                    set(handles.cluster0button,'value',newVal);
                    
                    %invokde the callback
                    matclust('cluster0button_Callback',handles.cluster0button,handles)
                else
                    %here is the handle for the radio button we want to toggle
                    radObj= figattrib.clustcontrol(clusNum,4);
                    
                    %toggle it
                    newVal = (get(radObj,'value')==0)+0;
                    set(radObj,'value',newVal);
                    
                    %invoke its callback
                    matclust('clustradiobutton2_Callback',radObj,handles);
                    
                    
                    
                end
            end
            
           
            
	end
end
%-----------------------------------------------------------
function highlight_ButtonDownFcn(hObject,handles)
%called if a hightlighted square is pressed
%because clicking on the square blocks the needed call for the
%clicked polygon, we need to route the buttondown_Fcn with the correct
%polygon handle.

global graphattrib;

hdata = get(hObject,'UserData');
polyindex = hdata.polyindex;
polyhandle = graphattrib.polyg(polyindex(1),polyindex(2),polyindex(3)).lines{polyindex(4)};
polygon_ButtonDownFcn(polyhandle,handles);

%-----------------------------------------------------------
function clearhighlight()
%un-highlights the highlighted cluster polygons

global figattrib;
global graphattrib;

graphattrib.squarehighlighton = 0;
highlightindex = [];
%set(figattrib.handles.mainPolygonMenu','Enable','off');
set(figattrib.handles.polydeletemenu,'Enable','off');
set(figattrib.handles.addPolygonFilterMainmenu,'Enable','off');
set(figattrib.handles.copyPolygonFilterMainmenu,'Enable','off');
try 
    set(graphattrib.currentpolyhighlight,'Marker','none'); %turn any previous highlighting off 
end

try
    highlightindex = get(graphattrib.currentpolyhighlight,'UserData');
end
if ~isempty(highlightindex)
    if ~iscell(highlightindex)
        try
            delete(graphattrib.polyg(highlightindex(1),highlightindex(2),highlightindex(3)).highlight{highlightindex(4)}); %clear any previous vertex highlight 
        end        
    else
        for i = 1:length(highlightindex)
            try
                delete(graphattrib.polyg(highlightindex{i}(1),highlightindex{i}(2),highlightindex{i}(3)).highlight{highlightindex{i}(4)}); %clear any previous vertex highlight 
            end
        end
    end
end

graphattrib.currentpolyhighlight = [];

%-----------------------------------------------------------------
function plotgraph(handles)
%plots the information on the desired axes (axes are controlled by the two listboxes)
%this is the main plotting function

global clustdata;
global figattrib;
global graphattrib;
global clustattrib;

timeclock = toc;

if ~isempty(clustdata.params)
    
	currx = get(handles.listbox1,'Value'); %the column number in the parameter matrix to plot on the x axis
	curry = get(handles.listbox2,'Value'); %the column number in the parameter matrix to plot on the y axis
	tempindex = [];
	numclustdata.params = size(clustdata.params,2);
	axes(graphattrib.graphwindow);
	s = get(graphattrib.graphwindow,'Position');
	
    set(graphattrib.graphwindow,'Color',graphattrib.backgroundcolor); %set the background color of the plot
	resolutionx = s(3)*graphattrib.resolutionfactor(1); %calculate the resolution of the plot image	
    resolutiony = s(4)*graphattrib.resolutionfactor(2);
	currviewbox = graphattrib.viewbox(currx,:,curry); %get the current view range of the desired data
	step = 4;
    
    %The first two colors in the colormap are the background color and the
    %color of cluster 0.  Then comes clusters 1 through ...
    colormap([graphattrib.backgroundcolor;clustattrib.cluster0attrib.color;figattrib.mixcolor]); 
                                                                                                            
	
	%assign cluster numbers for each cluster-- for coloring (1 is the background color, 2 is
	%the color of cluster 0, 3 is the color of the first cluster, ...)
	tempindex = [(1:length(clustdata.filteredpoints))' ones(size(clustdata.filteredpoints))];
	tempfilter = clustdata.filteredpoints; 
	graphattrib.nonhiddenpoints = true(length(clustdata.filteredpoints),1);
    graphattrib.blackedOutClusters = [];
    if (clustattrib.cluster0attrib.show) %assign color to cluster 0 if it is currently visible
        tempindex(:,2) = 2;
    end
	for i = clustattrib.clustersOn'    
        if ((clustattrib.hiddenclusters(i,1)==0)&(clustattrib.hiddenclusters(i,2)==0)) %don't hide cluster
            tempindex(clustattrib.clusters{i}.index,2) = i+2; 
        elseif (clustattrib.hiddenclusters(i,1)==1) %hide cluster by turning points black
            tempfilter(clustattrib.clusters{i}.index) = 0;
            graphattrib.nonhiddenpoints(clustattrib.clusters{i}.index) = false;
            graphattrib.blackedOutClusters = [graphattrib.blackedOutClusters i]; 
        end
        
        pointcount = length(find(tempfilter(clustattrib.clusters{i}.index)));
        set(figattrib.clustcontrol(i,2),'TooltipString',['Cluster ',num2str(i),': ',num2str(pointcount),' points']);
    end
  
    tempfilter(find(tempindex(:,2) == 1)) = 0;   
	tempindex = tempindex(find(tempfilter),:);  %remove the currently filtered points
	
    %if there is currently no data open
    if ~(clustattrib.nodata)
        pointcount = length(find(tempindex(:,2)==2));
        set(handles.cluster0button,'TooltipString',['Cluster 0: ', num2str(pointcount),' points']);
    else
        set(handles.cluster0button,'TooltipString','Cluster 0');
    end
          
    %orderindex is a mex program-- it reorders tempindex so that the 
    %2's are first (cluster 0), 1's next, then the colors are ordered to
    %the user's preference
    tempindex = orderindex(tempindex',clustattrib.clustersOn);
    tempindex = tempindex';
       
	%pick the two dimensions of data to plot
    if (currx > clustdata.origparams)
        tmprotdata = clustdata.params(:,clustdata.rotation(currx).params);
        tmprotdata(:,4) = 1;
        
        
        projection = (clustdata.rotation(currx).matrix*(tmprotdata'))';
        %plotspace = projection(tempindex(:,1),[2 1])./repmat(projection(tempindex(:,1),4),[1 2]);
        plotspace = projection(tempindex(:,1),[2 1]);
        
        clear tmpdata;
        clear projection;
    else
    
        plotspace = clustdata.params(tempindex(:,1),[curry currx]);
    end
     
    try    
        
        tempindex(end+1,1:2) = [0 2];  
        
        if (currx > clustdata.origparams)
            currentdatarange = clustdata.rotation(currx).datarange(1:2,[2 1]);
        else
            currentdatarange = [clustdata.datarange(1,curry) clustdata.datarange(1,currx);clustdata.datarange(2,curry) clustdata.datarange(2,currx)];
        end
        
        graphattrib.currentdatarange = currentdatarange;
        plotspace(end+1,1:2) = [currentdatarange(2,1) currentdatarange(2,2)];
         
        %make the minimum value fall on (1,1) because negative values are not
        %allowed when making an image
        
        plotspace(:,1) = plotspace(:,1) - (currentdatarange(1,1));
        plotspace(:,2) = plotspace(:,2) - (currentdatarange(1,2));

        
        %the data points are changed for the correct resolution

        %currviewbox is [xmin xmax ymin ymax]-- numbers between 0 and 1 describing 
        %the percentage of data to be shown on each axis.  The resolution of the
        %data is increased with smaller window sizes.

        %the data is also normalized by the full datarange in each axis.
        %summary of calculation ---- plotspace = round((1/window)*resolution*normalizeddata)

        warning off MATLAB:divideByZero;
        graphscale = diag([(1/(currviewbox(4)-currviewbox(3)))*resolutiony/((currentdatarange(2,1)-currentdatarange(1,1))) ; (1/(currviewbox(2)-currviewbox(1)))*resolutionx/((currentdatarange(2,2)-currentdatarange(1,2)))]);
        plotspace = round(plotspace*graphscale);
         
        if floor((1/graphattrib.resolutionfactor(1))/2) %for large-point mode, the image should be shifted slightly
            plotspace = plotspace-floor((1/graphattrib.resolutionfactor(1))/2);
        end
        plotspace(find(plotspace(:)<1)) = 1; %smallest possible point is now (1,1)
         
         
        
        %convert the list of points to a sparse grid and
        %get rid of duplicate points
        %this keeps the last occurance of each duplicate
        %because cluster0 is at the top, it becomes overwritten by duplicates
        %one of two methods is used, depending on the current view resolution (for
        %speed)
        
        %[(1/(currviewbox(4)-currviewbox(3)))*resolutiony (1/(currviewbox(2)-currviewbox(1)))*resolutionx]
        if ~isempty(plotspace)
            if (max([(1/(currviewbox(4)-currviewbox(3)))*resolutiony (1/(currviewbox(2)-currviewbox(1)))*resolutionx])<5000)                
                spplotspace = zeros(max(plotspace(:,1)),max(plotspace(:,2)));               
                indconv = ((plotspace(:,2)-1)*max(plotspace(:,1)))+plotspace(:,1);                
                spplotspace(indconv) = tempindex(:,2); %assign figattrib.colors to the points               
            else	
                [plotspace, I,J] = unique(plotspace,'rows');
                %assign colors to the points
                plotspace(:,3) = tempindex(I,2);
                %convert the list of points to a sparse grid
                try
                    spplotspace = spconvert(plotspace);
                catch
                    graphattrib.viewbox = graphattrib.oldviewbox;
                    return
                end
            end
        else
            spplotspace = sparse(resolutiony,resolutionx);  %show no points if the filters have taken out all points
        end
    
    catch
        spplotspace = sparse(resolutiony,resolutionx);  
    end
     
    %newaxis is a scaled version of currviewbox (percentage->pixel number) 
    newaxis = round([(1/(currviewbox(2)-currviewbox(1)))*resolutionx*[currviewbox(1) currviewbox(2)] (1/(currviewbox(4)-currviewbox(3)))*resolutiony*[currviewbox(3) currviewbox(4)]]);

    %viewrage is [theoreticalxresolution, theoreticalyresolution, viewwindow(4
    %values)]  
    viewrange = [(1/(currviewbox(2)-currviewbox(1)))*resolutionx (1/(currviewbox(4)-currviewbox(3)))*resolutiony max([newaxis(1) 1]) min([newaxis(2) size(spplotspace,2)]) max([newaxis(3) 1]) min([newaxis(4) size(spplotspace,1)])]; 


    xrange = viewrange(4)-viewrange(3);
    yrange = viewrange(6)-viewrange(5);
    if (viewrange(5) < yrange)
        yrange = viewrange(5)-1;
    end
    if (viewrange(3) < xrange)
        xrange = viewrange(3)-1;
    end


    figattrib.spplotspace = spplotspace;

    hold on
    if (graphattrib.Himage == 0) %no image object exists yet (during startup)
        graphattrib.Himage = image(spplotspace,'EraseMode','normal');
        set(graphattrib.graphwindow,'YDir','normal');
        hold on    
        axis([newaxis(1)-viewrange(3) newaxis(2)-viewrange(3) newaxis(3)-viewrange(5) newaxis(4)-viewrange(5)]);           
    else
        set(graphattrib.Himage,'CData',spplotspace(viewrange(5):viewrange(6),viewrange(3):viewrange(4)));
        axis([newaxis(1)-viewrange(3) newaxis(2)-viewrange(3) newaxis(3)-viewrange(5) newaxis(4)-viewrange(5)]);                  
    end
        
    set(graphattrib.graphwindow,'UserData',viewrange);

    %change the axis tick numbers to match the plotted parameters
    XTick = get(graphattrib.graphwindow,'XTick');
    YTick = get(graphattrib.graphwindow,'YTick');
    newXTick = [1 resolutionx/4  resolutionx/2 3*resolutionx/4 resolutionx];
    newYTick = [1 resolutiony/4  resolutiony/2 3*resolutiony/4 resolutiony];
    set(graphattrib.graphwindow,'XTick',newXTick);
    set(graphattrib.graphwindow,'YTick',newYTick);
    XTick = newXTick;
    YTick = newYTick;
    xrange = ((abs(diff(currentdatarange(:,2))))*[currviewbox(1);currviewbox(2)])+ (currentdatarange(1,2));
    yrange = ((abs(diff(currentdatarange(:,1))))*[currviewbox(3);currviewbox(4)])+ (currentdatarange(1,1));
    xpercentages = (XTick-min(XTick))/(max(XTick)-min(XTick));
    ypercentages = (YTick-min(YTick))/(max(YTick)-min(YTick));
    newxlabels = round((((abs(diff(xrange)))*xpercentages)+xrange(1))*10)/10;
    if (currx == 1)
        newxlabels = timetrans(newxlabels',clustdata.UnitsPerSec,1)';
    end
    newylabels = round((((abs(diff(yrange)))*ypercentages)+yrange(1))*10)/10;
    if (curry == 1)
        newylabels = timetrans(newylabels',clustdata.UnitsPerSec,1)';
    end

    set(graphattrib.graphwindow,'XTickLabel',newxlabels);
    set(graphattrib.graphwindow,'YTickLabel',newylabels);
    set(graphattrib.graphwindow,'FontSize',7);

    %show the cluster polygons in the current axes 
    [findex, filterline] = findfilterindex(handles);
    for fNum = 1:length(findex)
        try
            for i = 1:length(graphattrib.polyg(currx,curry,findex(fNum)).vertices)
                try
                    graphattrib.polyg(currx,curry,findex(fNum)).vertices{i} = [(graphattrib.polyg(currx,curry,findex(fNum)).relativevertices{i}(:,1)*viewrange(1))-viewrange(3) (graphattrib.polyg(currx,curry,findex(fNum)).relativevertices{i}(:,2)*viewrange(2))-viewrange(5)];
                    set(graphattrib.polyg(currx,curry,findex(fNum)).lines{i},'XData',graphattrib.polyg(currx,curry,findex(fNum)).vertices{i}(:,1));
                    set(graphattrib.polyg(currx,curry,findex(fNum)).lines{i},'YData',graphattrib.polyg(currx,curry,findex(fNum)).vertices{i}(:,2));
                    set(graphattrib.polyg(currx,curry,findex(fNum)).lines{i},'Color',colorclip(figattrib.mixcolor(i,:)-[.2 .2 .2]));    
                    if (clustattrib.hiddenclusters(i,2))
                        set(graphattrib.polyg(currx,curry,findex(fNum)).lines{i},'Visible','off');    
                    else
                        set(graphattrib.polyg(currx,curry,findex(fNum)).lines{i},'Visible','on');
                    end
                end
            end
            if (graphattrib.plothighlightclear) %clear the polygon highlight?
                clearhighlight;
            end 
        end
    end
    graphattrib.plothighlightclear = 1;

else %the params matrix is empty-- no data to plot
    
    axes(graphattrib.graphwindow);
    spplotspace = sparse(100,100); 
    graphattrib.Himage = image(spplotspace,'EraseMode','normal');
    set(graphattrib.graphwindow,'YDir','normal');
    figattrib.spplotspace = spplotspace;
end
try
    if ~isempty(figattrib.rotationwindowHandle)
        rotationfunc('filter_Callback',0,guidata(figattrib.rotationwindowHandle))
    end
end
set(handles.statusText,'String','Ready');
set(handles.statusText,'UserData','Ready');
%-------------------------------------------------------------



%Functions for figure control
%----------------------------------------------
%-----------------------------------------------
function figure1_ResizeFcn(hObject,handles)

%called whenever the user resizes the main figure window

global graphattrib;
global figattrib;



sidePanelWidth = figattrib.sidePanelWidth; %the width of the side control panel
filterBoxDivider = figattrib.filterBoxDivider; %the division point between the two filter list boxes (range between 0 and 1) 

if (sidePanelWidth < 265)
    sidePanelWidth = 265;
    figattrib.sidePanelWidth = 265;
end
if ((filterBoxDivider < 0)|(filterBoxDivider > 1))
    filterBoxDivider = .5;
    figattrib.filterBoxDivider = .5;
end

releasepolydraw(handles);
figsize = get(hObject,'Position');
currdir = pwd;
cd(figattrib.foldername);
load matclust_defaults;
matclust_defaults.FigPosition = figsize;
save('matclust_defaults','matclust_defaults');
cd(currdir);

topyloc = (2*figsize(4)/3)-24;
buttonsep = 28;
set(handles.polygonbutton,'Position',[2+(buttonsep*0) figsize(4)-21 30 22]);
set(handles.squarebutton,'Position',[2+(buttonsep*1) figsize(4)-21 30 22]);
set(handles.arrowbutton,'Position',[2+(buttonsep*2) figsize(4)-21 30 22]);
set(handles.magbutton,'Position',[2+(buttonsep*3) figsize(4)-21 30 22]);
set(handles.handbutton,'Position',[2+(buttonsep*4) figsize(4)-21 30 22]);
set(handles.wandbutton,'Position',[2+(buttonsep*5) figsize(4)-21 30 22]);
set(handles.wandSlider,'Position',[7+(buttonsep*6) figsize(4)-22 200 20]);
set(handles.wandSliderLabel,'Position',[215+(buttonsep*6) figsize(4)-20 200 15]);

set(handles.zoominbutton,'Position',[figsize(3)-sidePanelWidth-33 figsize(4)-20 30 22]);
set(handles.zoomoutbutton,'Position',[figsize(3)-sidePanelWidth-63 figsize(4)-20 30 22]);

set(graphattrib.graphwindow,'Position',[40 35 figsize(3)-sidePanelWidth-45 figsize(4)-56]);
set(handles.frame2,'Position',[figsize(3)-sidePanelWidth+5 1*figsize(4)/3 sidePanelWidth figsize(4)/3]);

buttonlistlength = 11*length(figattrib.mixcolor);
yloc = topyloc-get(handles.clustslider,'Value')+buttonlistlength;
xloc = figsize(3)-sidePanelWidth+23;
oldloc = get(figattrib.clustcontrol(1,1),'Position');
oldloc = oldloc(1:2);
locdiff = [xloc yloc]-oldloc;


currloc = get(figattrib.clustcontrol,'Position');
for i = 1:length(currloc)
     currloc{i}(1:2) = currloc{i}(1:2)+locdiff;
end
set(figattrib.clustcontrol,{'Position'},currloc);

spacer = (sidePanelWidth-265)/2;
set(handles.cluster0button,'Position',[figsize(3)-sidePanelWidth+150 (2*figsize(4)/3)-23 sidePanelWidth-151 20]);
set(handles.hideallonbutton,'Position',[figsize(3)-sidePanelWidth+193+spacer 2*(figsize(4)/3)-45 25 17]);
set(handles.hidealloffbutton,'Position',[figsize(3)-sidePanelWidth+217+spacer 2*(figsize(4)/3)-45 45 17]);
set(handles.hidealltext,'Position',[figsize(3)-sidePanelWidth+148+spacer 2*(figsize(4)/3)-49 44 20]);
set(handles.excludeallonbutton,'Position',[figsize(3)-sidePanelWidth+193+spacer 2*(figsize(4)/3)-65 25 17]);
set(handles.excludealloffbutton,'Position',[figsize(3)-sidePanelWidth+217+spacer 2*(figsize(4)/3)-65 45 17]);
set(handles.excludealltext,'Position',[figsize(3)-sidePanelWidth+148+spacer 2*(figsize(4)/3)-69 44 20]);

set(handles.clusterinfo,'Position',[figsize(3)-sidePanelWidth+150 (figsize(4)/3)+3 sidePanelWidth-151 (figsize(4)/3)-75]);
set(handles.clustslider,'Position',[figsize(3)-sidePanelWidth+6 figsize(4)/3+1 17 figsize(4)/3-2]);
set(handles.frame1,'Position',[figsize(3)-sidePanelWidth+5 (2*figsize(4)/3)-1 sidePanelWidth (figsize(4)/3)+10]);
set(handles.frame3,'Position',[figsize(3)-sidePanelWidth+5 20 sidePanelWidth (figsize(4)/3)-18]);
set(handles.timefilterText,'Position',[figsize(3)-(sidePanelWidth*(1-filterBoxDivider))-80 (figsize(4)/3)-22 120 20]);
set(handles.otherfilterText,'Position',[figsize(3)-82 (figsize(4)/3)-22 120 20]);
set(handles.timefilterList,'Position',[figsize(3)-sidePanelWidth+8 25 (sidePanelWidth*filterBoxDivider)-7.5 (figsize(4)/3)-40]);
set(handles.otherfilterList,'Position',[figsize(3)-(sidePanelWidth*(1-filterBoxDivider))+6 25 (sidePanelWidth*(1-filterBoxDivider))-7.5 (figsize(4)/3)-40]);
set(handles.xText,'Position',[figsize(3)-(sidePanelWidth/2)-82 (figsize(4))-19 120 15]);
set(handles.yText,'Position',[figsize(3)-82 (figsize(4))-19 120 15]);
set(handles.listbox1,'Position',[figsize(3)-sidePanelWidth+8 (2*figsize(4)/3)+1 (sidePanelWidth/2)-7.5 (figsize(4)/3)-20]);
set(handles.listbox2,'Position',[figsize(3)-(sidePanelWidth/2)+6 (2*figsize(4)/3)+1 (sidePanelWidth/2)-7.5 (figsize(4)/3)-20]);
set(handles.statusBox,'Position',[0 0 figsize(3)+5 20]);
set(handles.statusText,'Position',[0 0 figsize(3) 18]);
set(handles.timeFiltOptButton,'Position',[figsize(3)-sidePanelWidth+10 (figsize(4)/3)-16 32 17]);

figattrib.lastredraw = toc;

if (strcmp(get(figattrib.timers,'Running'),'on'))
    stop(figattrib.timers);
end
start(figattrib.timers);


%------------------------------------------------------
function figure1_CloseRequestFcn(hObject,handles)

global figattrib;

allclosed = 0;
go = 0;

while (~allclosed)
    if ~isempty(figattrib.openfiles)
        name1 = figattrib.openfiles{figattrib.currentopenfile};
        closeMenu_Callback(hObject,handles);
        try
            if strcmp(name1,figattrib.openfiles{figattrib.currentopenfile})
                break
            end
        end
    else
        allclosed = 1;
        go = 1;
    end
end

if (go)
	currdir = pwd;
	cd(figattrib.datafoldername);
	if (ispc)
        system('rd /S /Q M_dataopen');
	elseif (isunix)
        system('rm -r -f M_dataopen');
	end
	%mkdir('M_dataopen');
	cd(currdir);
	try
        stop(figattrib.timers);
        delete(figattrib.timers);
    end
    
    delete(figattrib.handles.figure1);
    try
        if ~isempty(figattrib.rotationwindowHandle)
            delete(figattrib.rotationwindowHandle);
        end
    end
    clear global figattrib;
    clear global clustattrib;
    clear global clustdata;
    clear global graphattrib; 

end
%-------------------------------------------------------
function redrawTimerFcn()
%timer function that redraws the axes after a resize (there may be a better
%way to do this...)
global figattrib;
global clustattrib;

if (toc-figattrib.lastredraw)>.05
    plotgraph(figattrib.handles);
    stop(figattrib.timers);
end
%-------------------------------------------------------------
function findInPoints(a1,a2,findex,clustnum, clustexcludes)
%calulates all of the data points that are within the polygons for the
%specified cluster. This function can be called with 4 or 5 inputs.  Also,
%if the first three inputs are empty matrices, then then it is assumed that
%no polygon has been created or changed, and that only the cluster excludes
%need to be recalculated.

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;

set(figattrib.handles.statusText,'String','Finding points inside shape...');
set(figattrib.handles.statusText,'UserData','Finding points inside shape...');

if (nargin == 4)
    %filteredpoints = tmpfilter & tmpfilter2 & graphattrib.nonhiddenpoints; 
    %if any clusters are blacked out, exclude them too
    clustExcludePoints = false(size(clustdata.timefilters,1),1);
elseif (nargin == 5) %if called with clustexcludes
    clustExcludePoints = false(size(clustdata.timefilters,1),1);
    for i = clustexcludes
        clustExcludePoints(clustattrib.clusters{i}.index) = true;
    end
end

if (~isempty(a1) & ~isempty(a2) & ~isempty(findex)) %if the first three inputs are not empty, it means a new polygon was added in the current axes
    newpolygon = 1;
else
    newpolygon = 0;
end

if (newpolygon) %a polygon was created/changed, so calculate the new points for the cluster
    filts = clustattrib.filterindex{findex};
    tmpfilter = fastandbit(clustdata.timefilters,filts(1)); %temporal filters
    tmpfilter2 = fastandbit(clustdata.otherfilters,filts(2:end)); %other filters
    filteredpoints = tmpfilter & tmpfilter2;
    if (nargin == 4)
        %filteredpoints = tmpfilter & tmpfilter2 & graphattrib.nonhiddenpoints; %if any clusters are blacked out, exclude them too
        clustattrib.dependencies(clustnum, 1:end) = false;
        clustattrib.dependencies(clustnum, graphattrib.blackedOutClusters) = true;
    elseif (nargin == 5) %if called with clustexcludes
        %filteredpoints = tmpfilter & tmpfilter2 & graphattrib.nonhiddenpoints & (~clustExcludePoints); %do not include these additional clusters
        clustattrib.dependencies(clustnum, 1:end) = false;
        clustattrib.dependencies(clustnum, clustexcludes) = true;
        clustattrib.dependencies(clustnum, graphattrib.blackedOutClusters) = true;
    end

    in = find(filteredpoints);
    %in = clustattrib.cluster0;  %all points allowed in current filter(s)

    if ~isempty(clustattrib.takenpolys) %if the current cluster had a box for the current axes and filters, that polygon is removed from takenpolys
        try

            removepoly = clustattrib.clusters{clustnum}.polyindex( ...
                find((clustattrib.clusters{clustnum}.polyindex(:,1)==a1)&(clustattrib.clusters{clustnum}.polyindex(:,2)==a2)&(clustattrib.clusters{clustnum}.polyindex(:,3)==findex)),4);
            clustattrib.takenpolys = sort(setdiff(clustattrib.takenpolys,removepoly));

        end
    end
    if ~isempty(clustattrib.takenpolys)  %the number of this polygon is the first available number in takenpolys
        polynum = min(setdiff(1:max(clustattrib.takenpolys)+1,clustattrib.takenpolys));
        clustattrib.takenpolys = sort([clustattrib.takenpolys polynum]);
    else
        polynum = 1;
        clustattrib.takenpolys = 1;
    end
    memvector = ceil(polynum/32);
    if (size(clustattrib.clusters{clustnum}.defineaxes,1) == 0)
        in = [];
    end
    
    if (a1 > clustdata.origparams)
        tmprotdata = clustdata.params(:,clustdata.rotation(a1).params);
        tmprotdata(:,4) = 1;
        
        
        projection = (clustdata.rotation(a1).matrix*(tmprotdata'))';
        plotspace = projection(:,[1 2]);
        
        clear tmpdata;
        clear projection;
    else
    
        plotspace = clustdata.params(:,[a1 a2]);
    end
    
    if (a1 > clustdata.origparams)
        currentdatarange = clustdata.rotation(a1).datarange(1:2,[2 1]);
    else
        currentdatarange = [clustdata.datarange(1,a2) clustdata.datarange(1,a1);clustdata.datarange(2,a2) clustdata.datarange(2,a1)];
    end
        
    xpoints = ((plotspace(:,1)-currentdatarange(1,2)))/(currentdatarange(2,2)-currentdatarange(1,2));
    ypoints = ((plotspace(:,2)-currentdatarange(1,1)))/(currentdatarange(2,1)-currentdatarange(1,1));
    
    %xpoints = ((plotspace(:,1)-graphattrib.currentdatarange(1,2)))/(graphattrib.currentdatarange(2,2)-graphattrib.currentdatarange(1,2));
    %ypoints = ((plotspace(:,2)-graphattrib.currentdatarange(1,1)))/(graphattrib.currentdatarange(2,1)-graphattrib.currentdatarange(1,1));

    %xpoints = ((clustdata.params(:,a1)-clustdata.datarange(1,a1)))/(clustdata.datarange(2,a1)-clustdata.datarange(1,a1));
    %ypoints = ((clustdata.params(:,a2)-clustdata.datarange(1,a2)))/(clustdata.datarange(2,a2)-clustdata.datarange(1,a2));

    vertices = graphattrib.polyg(a1,a2,findex).relativevertices{clustnum};
    xvert = vertices(:,1);
    yvert = vertices(:,2);
    excluding = false(length(clustdata.params(:,1)),1);
    including = false(length(clustdata.params(:,1)),1); %default assumption: all points are outside the polygon
    excluding(in) = true;  %default assumption: all points in the current filter are excluded by this polygon (note: points not in the filter are not excluded)
    in = in(find(fastinpoly(xpoints(in),ypoints(in),xvert,yvert))); %find the points of the current filter that are inside the polygon

    excluding(in) = false;  %points within polygon are not excluded
    including(in) = true; %points within polygon are included
    if (size(clustattrib.pointexclude,2)<memvector) %we have to add a new column to pointexclude and pointinclude for every 32 polygons
        clustattrib.pointexclude(:,end+1) = 0;
        clustattrib.pointinclude(:,end+1) = 0;
    end
    paramlength = length(clustdata.params(:,1));

    clustattrib.pointexclude(1:paramlength,memvector) = fastbitset(clustattrib.pointexclude(1:paramlength,memvector),polynum-(32*(memvector-1)),excluding);  %save the bit-wise info
    clustattrib.pointinclude(1:paramlength,memvector) = fastbitset(clustattrib.pointinclude(1:paramlength,memvector),polynum-(32*(memvector-1)),including);

    try
        clustattrib.clusters{clustnum}.polyindex = clustattrib.clusters{clustnum}.polyindex( ...
            find((clustattrib.clusters{clustnum}.polyindex(:,1)~=a1)|(clustattrib.clusters{clustnum}.polyindex(:,2)~=a2) ...
            |(clustattrib.clusters{clustnum}.polyindex(:,3)~=findex)),:);
        clustattrib.clusters{clustnum}.polyindex = [clustattrib.clusters{clustnum}.polyindex;[a1 a2 findex polynum]];
    catch
        clustattrib.clusters{clustnum}.polyindex = [a1 a2 findex polynum];
    end
end

in = fastorbit(clustattrib.pointexclude,clustattrib.clusters{clustnum}.polyindex(:,4)); %which points are excluded at any polygon
in2 = fastorbit(clustattrib.pointinclude,clustattrib.clusters{clustnum}.polyindex(:,4)); %which points have been included at any polygon
in3 = false(length(clustdata.params(:,1)),1); %which points are excluded individually

try
	if ~isempty(clustattrib.eventeditindex{clustnum})
        in3(clustattrib.eventeditindex{clustnum}(:,1)) = true; %if any points were excluded individually, it will be stored in clustattrib.eventeditindex
	end
end

if (newpolygon)
    clustattrib.clusters{clustnum}.index = uint32(find(~in & in2 & ~in3 & graphattrib.nonhiddenpoints & ~clustExcludePoints));
else
    %if the polygons for the cluster are unchanged, only the cluster
    %excludes need to be calculated
    clustattrib.clusters{clustnum}.index = uint32(find(~in & in2 & ~in3 & ~clustExcludePoints));
end

%any clusters that are dependent on this cluster must also be redefined
%(this will work its way down the tree of dependencies)
for i = find(clustattrib.dependencies(:,clustnum))'
    clustexcludes = find(clustattrib.dependencies(i,:));
    findInPoints([],[],[],i, clustexcludes);    
end
%------------------------------------------------------------
function psizeMenuFcn(hObject,handles)
%changes the size of the spike points

global graphattrib;

choice = get(hObject,'UserData');

if (choice == 1)
    graphattrib.resolutionfactor = [1 1];
    set(handles.smallpsizeMenu,'Checked','on');
    set(handles.largepsizeMenu,'Checked','off');
elseif (choice == 2)
    graphattrib.resolutionfactor = [.5 .5];
    set(handles.smallpsizeMenu,'Checked','off');
    set(handles.largepsizeMenu,'Checked','on');
end

plotgraph(handles)
addnewstate('change point size',handles); 
%------------------------------------------------------------
function tofrontFcn(hObject,handles)
%callback function to change the display order of the cluster
global clustattrib;

choice = get(hObject,'UserData');
notindex = find(clustattrib.clustersOn ~= choice);
isindex = find(clustattrib.clustersOn == choice);

if ((~isempty(notindex))&(~isempty(isindex)))
    clustattrib.clustersOn = [clustattrib.clustersOn(notindex);clustattrib.clustersOn(isindex)];
end
plotgraph(handles)
addnewstate('change order',handles); 
%------------------------------------------------------------
function tobackFcn(hObject,handles)
%callback function to change the display order of the cluster
global clustattrib;

choice = get(hObject,'UserData');
notindex = find(clustattrib.clustersOn ~= choice);
isindex = find(clustattrib.clustersOn == choice);

if ((~isempty(notindex))&(~isempty(isindex)))
    clustattrib.clustersOn = [clustattrib.clustersOn(isindex);clustattrib.clustersOn(notindex)];
end
plotgraph(handles)
addnewstate('change order',handles); 
%------------------------------------------------------------
function forewardoneFcn(hObject,handles)
%callback function to change the display order of the cluster
global clustattrib;

choice = get(hObject,'UserData');
notindex = find(clustattrib.clustersOn ~= choice);
isindex = find(clustattrib.clustersOn == choice);

if ((~isempty(notindex))&(~isempty(isindex)))
    if (isindex < length(clustattrib.clustersOn))
        clustattrib.clustersOn(isindex:isindex+1) = [clustattrib.clustersOn(isindex+1);clustattrib.clustersOn(isindex)];
    end
end
plotgraph(handles)
addnewstate('change order',handles); 
%------------------------------------------------------------
function backoneFcn(hObject,handles)
%callback function to change the display order of the cluster
global clustattrib;

choice = get(hObject,'UserData');
notindex = find(clustattrib.clustersOn ~= choice);
isindex = find(clustattrib.clustersOn == choice);

if ((~isempty(notindex))&(~isempty(isindex)))
    if (isindex > 1)
        clustattrib.clustersOn(isindex-1:isindex) = [clustattrib.clustersOn(isindex);clustattrib.clustersOn(isindex-1)];
    end
end
plotgraph(handles)
addnewstate('change order',handles); 
%------------------------------------------------------------

function deletepoly(a1,a2,findex,clustnum)
%called to delete a cluster box in the current axes

global graphattrib;
global clustattrib;
global figattrib;
global clustdata;

try 
    delete(graphattrib.polyg(a1,a2,findex).lines{clustnum});
    graphattrib.polyg(a1,a2,findex).lines{clustnum} = [];
end
%turn highlighting off
try 
    delete(graphattrib.polyg(a1,a2,findex).highlight{clustnum});             
end
%remove the record of a box in these axes if one exists
try
    clustattrib.clusters{clustnum}.defineaxes = setdiff(clustattrib.clusters{clustnum}.defineaxes,[a1 a2 findex],'rows');
catch
    clustattrib.clusters{clustnum}.defineaxes = [];
end
graphattrib.polyg(a1,a2,findex).vertices{clustnum} = [];
graphattrib.polyg(a1,a2,findex).relativevertices{clustnum} = [];
clustattrib.clusters{clustnum}.box{a1,a2,findex} = [];

%if there was only one cluster box for this cluster, then turn cluster off
%and remove any single exclude points

if (isempty(clustattrib.clusters{clustnum}.defineaxes))
    turnoffcluster(clustnum,figattrib.handles);
    clustattrib.eventeditindex{clustnum} = [];
end

try   
    removepoly = clustattrib.clusters{clustnum}.polyindex( ... 
        find((clustattrib.clusters{clustnum}.polyindex(:,1)==a1)&(clustattrib.clusters{clustnum}.polyindex(:,2)==a2)&(clustattrib.clusters{clustnum}.polyindex(:,3)==findex)),4);
    clustattrib.takenpolys = sort(setdiff(clustattrib.takenpolys,removepoly));
    clustattrib.clusters{clustnum}.polyindex = clustattrib.clusters{clustnum}.polyindex( ... 
        find((clustattrib.clusters{clustnum}.polyindex(:,1)~=a1)|(clustattrib.clusters{clustnum}.polyindex(:,2)~=a2)|(clustattrib.clusters{clustnum}.polyindex(:,3)~=findex)) ...
         ,:);
       
    clustexcludes = find(clustattrib.dependencies(clustnum,:));
    %we need to recalculate the points for the cluster
    findInPoints([],[],[],clustnum, clustexcludes);    
                                                                               
catch
    clustattrib.clusters{clustnum}.polyindex = [];
    clustattrib.clusters{clustnum}.index = []
end
%-------------------------------------------------------
function FilterPoints()
%calculates which points pass the currently activated time filters and
%'other' filters

global clustdata;

tmpfilter = fastorbit(clustdata.timefilters,find(clustdata.timefiltersOn));
tmpfilter2 = fastandbit(clustdata.otherfilters,find(clustdata.otherfiltersOn));

clustdata.filteredpoints = tmpfilter & tmpfilter2;
%---------------------------------------------------------
function [linepick, filterline] = findfilterindex(handles)

global clustattrib;
global clustdata;

t = clustdata.timefiltersOn(:)';
t = find(t);

o = clustdata.otherfiltersOn(:)';
o = find(o);

filterline = [];
linepick = [];
for timefilts = 1:length(t)
	filterline = [filterline; [t(timefilts) o]];
	linepick(timefilts) = -1;
	if ~isempty(clustattrib.filterindex)
        for i = 1:length(clustattrib.filterindex)
            if (isequal(clustattrib.filterindex{i},filterline(timefilts,:)))
                linepick(timefilts) = i;
                break;
            end
        end
	end
end
%--------------------------------------------------------
function linepick = addfilterindex(filterline)

global clustattrib;

linepick = -1;
if ~isempty(clustattrib.filterindex)
    for i = 1:length(clustattrib.filterindex)
        if isempty(clustattrib.filterindex{i})
            linepick =  i;
            break;
        end
    end
end
if (linepick == -1)
    linepick = length(clustattrib.filterindex)+1;
end
clustattrib.filterindex{linepick} = filterline;
%--------------------------------------------------------
function timefilterdelete(hObject,handles)

%called if the user deletes a time filter
%this also deletes all polygons using that filter

global clustdata;
global graphattrib;
global clustattrib;

reply = questdlg('Warning: this will delete all cluster boxes using this filter.','Warning','Cancel','Delete filter','Delete filter');
if strcmp(reply, 'Delete filter')
	releasepolydraw(handles);
    names = get(handles.timefilterList,'String');
	currval = get(handles.timefilterList,'Value');
	
	names{currval} = [num2str(currval),' '];   
	clustdata.timefilterranges(currval,1:2) = [0 0];
	memoryindex = find(clustdata.timefiltermemmap == currval);
	
    hidecurrentpolygons(handles); %hide currently displayed polygons   		
    clustdata.timefiltersOn(find(clustdata.otherfiltermemmap == currval)) = 0; %if the deleted filter was on, turn it off
    
    % if that was the only time filter on, then turn on alltimes filter
    if isempty(find(clustdata.timefiltersOn))
        tempname = names{1};
        tempname(min(strfind(tempname,' '))) = 187;
        names{1} = tempname;
        clustdata.timefiltersOn(1) = 1;
    end
	
	set(handles.timefilterList,'String',names);
    clustdata.timefilternames = names;
	%set(handles.timefilterList,'Value',newval);
	FilterPoints;
    showcurrentpolygons(handles);  %disply the polygons in the new filter
	
	%find which filter indeces used the deleted time filter
	linepick = [];
	if ~isempty(clustattrib.filterindex)
        for i = 1:length(clustattrib.filterindex)
            
            if (clustattrib.filterindex{i}(1) == memoryindex)
                linepick =  [linepick i];
                clustattrib.filterindex{i} = [];    
            end
        end
	end
	
	%delete the polygons using the filter index
	if ~isempty(linepick)
        for p = linepick
            for c = clustattrib.clustersOn'
                dpoly = min(find(clustattrib.clusters{c}.defineaxes(:,3) == p));
                while ~isempty(dpoly)
                    deletepoly(clustattrib.clusters{c}.defineaxes(dpoly,1),clustattrib.clusters{c}.defineaxes(dpoly,2),p,c);
                    dpoly = min(find(clustattrib.clusters{c}.defineaxes(:,3) == p));    
                end
            end
        end
	end
	clustdata.timefiltermemmap(memoryindex) = 0; % erase the index to the filter
	updateclustinfo(handles);
	cleartimefilter;
	plotgraph(handles);
    addnewstate('delete time filter',handles); 
end
%---------------------------------------------------------
function otherfilterdelete(hObject,handles)

%called if the user deletes a time filter
%this also deletes all polygons using that filter

global clustdata;
global graphattrib;
global clustattrib;

reply = questdlg('Warning: this will delete all cluster boxes using this filter.','Warning','Cancel','Delete filter','Delete filter');
if strcmp(reply, 'Delete filter')
	releasepolydraw(handles);
    names = get(handles.otherfilterList,'String');
	currval = get(handles.otherfilterList,'Value');
	
	names{currval} = [num2str(currval),' '];   	
	memoryindex = find(clustdata.otherfiltermemmap == currval);
    % 2012/01/24 Add by HL. Delete clustdata.otherfilterData when delete filter
    clustdata.otherfilterData{currval} = [];
    % ----- end-------------
			
    hidecurrentpolygons(handles); %hide currently displayed polygons		
    clustdata.otherfiltersOn(find(clustdata.otherfiltermemmap == currval)) = 0; %if the deleted filter was on, turn it off
	set(handles.otherfilterList,'String',names);
    clustdata.otherfilternames = names;
	%set(handles.otherfilterList,'Value',newval);
	FilterPoints;
	showcurrentpolygons(handles);  %disply the polygons in the new filter
    
	%find which filter indeces used the deleted filter
	linepick = [];
	if ~isempty(clustattrib.filterindex)
        for i = 1:length(clustattrib.filterindex)
            filtfind = find(clustattrib.filterindex{i}(2:end) == memoryindex);
            
            if (~isempty(filtfind))
                linepick =  [linepick i];
                clustattrib.filterindex{i} = [];    
            end
        end
	end
	
	%delete the polygons using the filter index
	if ~isempty(linepick)
        for p = linepick
            for c = clustattrib.clustersOn'
                dpoly = min(find(clustattrib.clusters{c}.defineaxes(:,3) == p));
                while ~isempty(dpoly)
                    deletepoly(clustattrib.clusters{c}.defineaxes(dpoly,1),clustattrib.clusters{c}.defineaxes(dpoly,2),p,c);
                    dpoly = min(find(clustattrib.clusters{c}.defineaxes(:,3) == p));    
                end
            end
        end
	end
	clustdata.otherfiltermemmap(memoryindex) = 0; % erase the index to the filter
    
    % 2012/01/24 Add by HL. Delete clustdata.otherfilterData when delete filter
    clustdata.filtermemmap = [clustdata.timefiltermemmap; clustdata.otherfiltermemmap];
    % ----- end-------------

	updateclustinfo(handles);
	clearotherfilter;
	plotgraph(handles);
    addnewstate('delete filter',handles); 
end
%--------------------------------------------------------
function resetfilters(resetval)

% resetval
% 0 -- reset all filter text to blank
% 1 -- keep text, update arrows
% 2 -- keep text, update arrows, if any filters on do not exist, update
% filtersOn values
global clustdata;
global graphattrib;
global clustattrib;
global figattrib;

if (resetval == 0)
    
    tmp = zeros(32,1);
    tmp(1) = 1;
    otherfiltersOn = tmp'; %for the 'other' filters, more than one filter can be on at the same time 
    timefiltersOn = tmp'; %which time filters are on
elseif (resetval == 1)
    timefiltersOn = clustdata.timefiltersOn;
    otherfiltersOn = clustdata.otherfiltersOn;
elseif (resetval == 2)
    timefiltersOn = clustdata.timefiltersOn;
    otherfiltersOn = clustdata.otherfiltersOn;
    
end

handles = figattrib.handles;
names = get(handles.timefilterList,'String');
releasepolydraw(handles);
timearrow = clustdata.timefiltermemmap(find(timefiltersOn)); 

if  ((timearrow ~= 0)|(resetval~=2))
    for i = 1:length(names)
        tempname = names{i};
        tempname(strfind(names{i},187)) = ' ';
        
        if ismember(i,timearrow)
            tempname(min(strfind(tempname,' '))) = 187; %add double arrows to highlight selection
        end
        names{i} = tempname;
        
    end
    clustdata.timefiltersOn = timefiltersOn; %finds the index to the filter information
else
    for i = 1:length(names)
        tempname = names{i};
        tempname(strfind(names{i},187)) = ' ';
        
        if (i == 1)
            tempname(min(strfind(tempname,' '))) = 187; %add double arrows to highlight selection
        end
        names{i} = tempname;
    end
    tmp = zeros(32,1);
    tmp(1) = 1;
    timefiltersOn = tmp'; 
 end
set(handles.timefilterList,'String',names);
clustdata.timefilternames = names;
names = get(handles.otherfilterList,'String');
otherarrows = clustdata.otherfiltermemmap(find(otherfiltersOn)); 

for i = 1:length(names)
    tempname = names{i};
    tempname(strfind(names{i},187)) = ' ';

    if ismember(i,otherarrows)
        tempname(min(strfind(tempname,' '))) = 187; %add double arrows to highlight selection
    end
    names{i} = tempname;
end        
set(handles.otherfilterList,'String',names);
if (resetval == 2)
    for i = 1:length(otherfiltersOn)
        if (otherfiltersOn(i)==1)
            if (clustdata.otherfiltermemmap(i) == 0)
                otherfiltersOn(i) = 0;
            end
        end
    end
end

clustdata.otherfilternames = names;

clustdata.otherfiltersOn = otherfiltersOn;          
%--------------------------------------------------------
function updateclustinfo(handles)

global clustattrib;

liststring = [];

a1 = get(handles.listbox1,'Value');
a2 = get(handles.listbox2,'Value');
[linepick, filterline] = findfilterindex(handles);

currentview = [a1 a2 linepick];
viewmatch = [];
try
	if ~isempty(clustattrib.clusters{clustattrib.currclust}.defineaxes)
        for i = 1:size(clustattrib.clusters{clustattrib.currclust}.defineaxes,1)
            tmpline = clustattrib.clusters{clustattrib.currclust}.defineaxes(i,:);
            if isequal(currentview,tmpline)
                viewmatch = i;
            end
            filters = clustattrib.filterindex{tmpline(3)};
            liststring{i} = ['x(',num2str(tmpline(1)),'),y(',num2str(tmpline(2)),'),t(',num2str(filters(1)),'),o(',num2str(filters(2:end)),')'];
        end    
    end
end
set(handles.clusterinfo,'Value',viewmatch);
set(handles.clusterinfo,'String',liststring);
if (ismember(clustattrib.currclust,clustattrib.clustersOn))
    set(handles.mainClustMenu,'Enable','on');
    for menuNum = 1:length(handles.mainClusterMenu)
        set(handles.mainClusterMenu(menuNum),'UserData',clustattrib.currclust)
    end
else
    set(handles.mainClustMenu,'Enable','off');
end
%----------------------------------------------------
function cleartimefilter()

global figattrib;

handles = figattrib.handles;
set(handles.timefilterList,'Value',[]);
set(handles.timeaddContext,'Enable','off');
set(handles.timedeleteContext,'Enable','off');
set(handles.timeeditContext,'Enable','off');
%----------------------------------------------------
function clearotherfilter()

global figattrib;

handles = figattrib.handles;
set(handles.otherfilterList,'Value',[]);
set(handles.otheraddContext,'Enable','off');
set(handles.otherdeleteContext,'Enable','off');
%-----------------------------------------------------
function clustdelete(clustnum, handles, draw)

global clustattrib;

releasepolydraw(handles);
def = clustattrib.clusters{clustnum}.defineaxes;
if ~isempty(def)
    for i = 1:size(def,1)

        deletepoly(def(i,1),def(i,2),def(i,3),clustnum);
    end
end
updateclustinfo(handles);    
if (draw)
    plotgraph(handles);
    addnewstate('delete cluster',handles); 
end
%-------------------------------------------------------
function clustfiltdelete(clustnum, handles)
%deletes all polygons for the current time filter(s)

global clustattrib;
global clustdata;

releasepolydraw(handles);
currval = get(handles.timefilterList,'Value');
memoryindex = find(clustdata.timefiltersOn); 
def = clustattrib.clusters{clustnum}.defineaxes;

%find which filter indeces used the deleted time filter
linepick = [];
if ~isempty(clustattrib.filterindex)
    for i = 1:length(clustattrib.filterindex)
        
        if ismember(clustattrib.filterindex{i}(1),memoryindex)
            linepick =  [linepick i];
            
        end
    end
end

%delete the polygons using the filter index
if (~isempty(def))&(~isempty(linepick))
    for i = 1:size(def,1)
        if ismember(def(i,3),linepick)
            deletepoly(def(i,1),def(i,2),def(i,3),clustnum);
        end    
    end
end
updateclustinfo(handles);
plotgraph(handles);
addnewstate('delete cluster',handles); 
%----------------------------------------------------
function clustwindelete(clustnum, handles)
%deletes the current polygon(s) for the cluster in the current axes

global clustattrib;

releasepolydraw(handles);
a1 = get(handles.listbox1,'Value');
a2 = get(handles.listbox2,'Value');
[linepick, filterline] = findfilterindex(handles);
def = clustattrib.clusters{clustnum}.defineaxes;
if ~isempty(def)   
    %if multiple time filters are selected, then it goes through and
    %deletes the box in each time filter
    for i = 1:length(linepick)
        deletepoly(a1,a2,linepick(i),clustnum);
    end
end
updateclustinfo(handles);
plotgraph(handles);
addnewstate('delete polygon(s)',handles); 
%-------------------------------------------------
function spexcludedelete(clustnum, handles)
%removes the 'exclude' status of all points in the cluster

global clustattrib;

clustattrib.eventeditindex{clustnum} = [];
in = fastorbit(clustattrib.pointexclude,clustattrib.clusters{clustnum}.polyindex(:,4)); %which points are excluded at any polygon
in2 = fastorbit(clustattrib.pointinclude,clustattrib.clusters{clustnum}.polyindex(:,4)); %which points have been included at any polygon
    
clustattrib.clusters{clustnum}.index = uint32(find(~in & in2));
plotgraph(handles);
addnewstate('delete single excludes',handles); 
%-----------------------------------------------
function copyclustMenuFcn(clustnum,handles)
%callback for the 'copy cluster' menu

global clustattrib
global figattrib

a = inputdlg('Enter an empty cluter number to copy into.','Copy cluster');
if isempty(a)
    return
end

a = str2num(a{1});
if (length(a)==1)&(~ismember(a,clustattrib.clustersOn))&(a < length(figattrib.mixcolor))
    
    clustcopy(clustnum,a,handles);
else
    errordlg('Copy error: enter a number for an unoccupied cluster ID');
end
    
%-----------------------------------------------
function copyepochMenuFcn(clustnum,handles)
%callback for the 'copy all polygons in time filter' menu

global clustattrib
global figattrib
global clustdata



b = inputdlg('Enter time filter to copy from.','Copy polygons');
if (~isempty(b))
    a = inputdlg('Enter time filter to copy into (can be a list).','Copy polygons');
    if isempty(a)
        return
    end
else
    return
end

names = get(handles.timefilterList,'String');
a = str2num(a{1});
a = round(a);
b = str2num(b{1});
b = round(b);

for i = 1:numel(a)
    
    if ((length(a(i))==1)&(~isempty(find(clustdata.timefiltermemmap == a(i))))&(a(i)>0)&(length(b)==1)&(~isempty(find(clustdata.timefiltermemmap == b)))&(b>0))
        
        epochcopy(clustnum,b,a(i),handles);
    else
        errordlg('Copy error: enter nonempty time filter IDs');
    end
end


% if ((length(a)==1)&(~isempty(find(clustdata.timefiltermemmap == a)))&(a>0)&(length(b)==1)&(~isempty(find(clustdata.timefiltermemmap == b)))&(b>0))
%     
%     epochcopy(clustnum,b,a,handles);
% else
%     errordlg('Copy error: enter nonempty time filter IDs'); 
% end


%---------------------------------------------
function clustcopy(clustnum,dest,handles)
%function to copy one cluster to another empty cluster number

global clustattrib;
global graphattrib;
global figattrib;
global clustdata;


def = clustattrib.clusters{clustnum}.defineaxes; 
releasepolydraw(handles);
if ~isempty(def)
    clustattrib.clusters{dest}.defineaxes = def;
    turnoncluster(dest,handles); %turn on cluster
    %set(figattrib.clustcontrol(dest,3:end),'Visible','on'); %turn on cluster control    
	if ~sum(ismember(clustattrib.clustersOn,dest))    
         clustattrib.clustersOn = [clustattrib.clustersOn;dest]; %add to 'on' list
    end
    %copy the cluster's dependencies
    clustattrib.dependencies(dest,:) = clustattrib.dependencies(clustnum,:);
    clustattrib.dependencies(:,dest) = clustattrib.dependencies(:,clustnum);
    for i = 1:size(def,1)
        
        a1 = def(i,1);
        a2 = def(i,2);
        findex = def(i,3);
        graphattrib.polyg(a1,a2,findex).vertices{dest} = graphattrib.polyg(a1,a2,findex).vertices{clustnum};
        graphattrib.polyg(a1,a2,findex).relativevertices{dest} = graphattrib.polyg(a1,a2,findex).relativevertices{clustnum};
        graphattrib.polyg(a1,a2,findex).lines{dest} = line(graphattrib.polyg(a1,a2,findex).vertices{clustnum}(:,1),graphattrib.polyg(a1,a2,findex).vertices{clustnum}(:,2));
        set(graphattrib.polyg(a1,a2,findex).lines{dest},'ButtonDownFcn','matclust(''polygon_ButtonDownFcn'',gcbo,guidata(gcbo))');
        
        set(graphattrib.polyg(a1,a2,findex).lines{dest},'UserData',[a1 a2 findex dest]);
        set(graphattrib.polyg(a1,a2,findex).lines{dest},'Color',colorclip(figattrib.mixcolor(dest,:)-[.2 .2 .2]));
        set(graphattrib.polyg(a1,a2,findex).lines{dest},'UIContextMenu',handles.polycontext);
        graphattrib.polyg(a1,a2,findex).type{dest} = 1;
        clustattrib.clusters{dest}.box{a1,a2} = graphattrib.polyg(a1,a2,findex).relativevertices{clustnum};    
        
        if ~isempty(clustattrib.takenpolys)  %the number of this polygon is the first available number in takenpolys
            polynum = min(setdiff(1:max(clustattrib.takenpolys)+1,clustattrib.takenpolys));
            clustattrib.takenpolys = sort([clustattrib.takenpolys polynum]);
        else
            polynum = 1;
            clustattrib.takenpolys = 1;
        end
        
        memvector = ceil(polynum/32);  
        
        clustattrib.clusters{dest}.polyindex(i,1:3) = clustattrib.clusters{clustnum}.polyindex(i,1:3);
        clustattrib.clusters{dest}.polyindex(i,4) = polynum;
        memvector2 = ceil(clustattrib.clusters{dest}.polyindex(i,4)/32);
       
        including = fastandbit(clustattrib.pointinclude,clustattrib.clusters{clustnum}.polyindex(i,4));
        excluding = fastandbit(clustattrib.pointexclude,clustattrib.clusters{clustnum}.polyindex(i,4));        
        
        if (size(clustattrib.pointexclude,2)<memvector) %we have to add a new column to pointexclude and pointinclude for every 32 polygons
            clustattrib.pointexclude(:,end+1) = 0;
            clustattrib.pointinclude(:,end+1) = 0;
        end
        
        clustattrib.pointexclude(1:length(clustdata.params(:,1)),memvector) = fastbitset(clustattrib.pointexclude(1:length(clustdata.params(:,1)),memvector),polynum-(32*(memvector-1)),excluding);  %save the bit-wise info
		clustattrib.pointinclude(1:length(clustdata.params(:,1)),memvector) = fastbitset(clustattrib.pointinclude(1:length(clustdata.params(:,1)),memvector),polynum-(32*(memvector-1)),including);
		
        curra1 = get(handles.listbox1,'UserData');
		curra2 = get(handles.listbox2,'UserData');
		[currfindex, filterline] = findfilterindex(handles);
		%hide the currently displayed polygons
		
        if ~isequal([curra1 curra2 currfindex],[a1 a2 findex])
            set(graphattrib.polyg(a1,a2,findex).lines{dest},'Visible','off');
        end
      
    end
    try
        if ~isempty(clustattrib.eventeditindex{clustnum})
            clustattrib.eventeditindex{dest} = clustattrib.eventeditindex{clustnum}; %if any points were excluded individually, it will be stored in clustattrib.eventeditindex
        end
    end
    clustattrib.clusters{dest}.index = clustattrib.clusters{clustnum}.index;
    plotgraph(handles);
    addnewstate('copy cluster',handles); 
end
%-----------------------------------------------
function copypoly(a1,a2,findex,clustnum,dest)
% allows the user to copy a polygon for cluster clustnum from time
% filter source to time filter dest.

global clustattrib;
global graphattrib;
global figattrib;
global clustdata;


sourcefilt = clustattrib.filterindex{findex};
destfilt = [dest sourcefilt(2:end)];
linepick = -1;

for i = 1:length(clustattrib.filterindex)
    if (isequal(clustattrib.filterindex{i},destfilt))
        linepick = i;
        break;
    end
end
if (linepick == -1)  %destination filter does not exist
    
    for i = 1:length(clustattrib.filterindex)
        if isempty(clustattrib.filterindex{i})
            linepick =  i;
            break;
        end
    end
end
if (linepick == -1)
    linepick = length(clustattrib.filterindex)+1;
end

clustattrib.filterindex{linepick} = destfilt;

graphattrib.polyg(a1,a2,linepick).vertices{clustnum} = graphattrib.polyg(a1,a2,findex).vertices{clustnum};
graphattrib.polyg(a1,a2,linepick).relativevertices{clustnum} = graphattrib.polyg(a1,a2,findex).relativevertices{clustnum};
graphattrib.polyg(a1,a2,linepick).lines{clustnum} = line(graphattrib.polyg(a1,a2,findex).vertices{clustnum}(:,1),graphattrib.polyg(a1,a2,findex).vertices{clustnum}(:,2));

set(graphattrib.polyg(a1,a2,linepick).lines{clustnum},'ButtonDownFcn','matclust(''polygon_ButtonDownFcn'',gcbo,guidata(gcbo))');

set(graphattrib.polyg(a1,a2,linepick).lines{clustnum},'UserData',[a1 a2 linepick clustnum]);
set(graphattrib.polyg(a1,a2,linepick).lines{clustnum},'Color',figattrib.mixcolor(clustnum,:));
set(graphattrib.polyg(a1,a2,linepick).lines{clustnum},'UIContextMenu',figattrib.handles.polycontext);
graphattrib.polyg(a1,a2,linepick).type{clustnum} = 1;

try
    clustattrib.clusters{clustnum}.defineaxes = setdiff(clustattrib.clusters{clustnum}.defineaxes,[a1 a2 linepick],'rows');
    clustattrib.clusters{clustnum}.defineaxes(end+1,:) = [a1 a2 linepick];
catch
    clustattrib.clusters{clustnum}.defineaxes = [a1 a2 linepick];
end

findInPoints(a1,a2,linepick,clustnum);
%-----------------------------------------------
function linepick = copypolyUniversal(a1,a2,findex,clustnum,destClust,destfilt)
% allows the user to copy a polygon for cluster clustnum from time
% filter source to time filter dest.

global clustattrib;
global graphattrib;
global figattrib;
global clustdata;

handles = figattrib.handles;
turnoncluster(destClust,handles); %turn on cluster
%set(figattrib.clustcontrol(dest,3:end),'Visible','on'); %turn on cluster control    
if ~sum(ismember(clustattrib.clustersOn,destClust))    
   clustattrib.clustersOn = [clustattrib.clustersOn;destClust]; %add to 'on' list
end
sourcefilt = clustattrib.filterindex{findex};
%destfilt = [dest sourcefilt(2:end)];
linepick = -1;

for i = 1:length(clustattrib.filterindex)
    if (isequal(clustattrib.filterindex{i},destfilt))
        linepick = i;
        break;
    end
end
if (linepick == -1)  %destination filter does not exist
    
    for i = 1:length(clustattrib.filterindex)
        if isempty(clustattrib.filterindex{i})
            linepick =  i;
            break;
        end
    end
end
if (linepick == -1)
    linepick = length(clustattrib.filterindex)+1;
end

clustattrib.filterindex{linepick} = destfilt;

graphattrib.polyg(a1,a2,linepick).vertices{destClust} = graphattrib.polyg(a1,a2,findex).vertices{clustnum};
graphattrib.polyg(a1,a2,linepick).relativevertices{destClust} = graphattrib.polyg(a1,a2,findex).relativevertices{clustnum};
graphattrib.polyg(a1,a2,linepick).lines{destClust} = line(graphattrib.polyg(a1,a2,findex).vertices{clustnum}(:,1),graphattrib.polyg(a1,a2,findex).vertices{clustnum}(:,2));

set(graphattrib.polyg(a1,a2,linepick).lines{destClust},'ButtonDownFcn','matclust(''polygon_ButtonDownFcn'',gcbo,guidata(gcbo))');

set(graphattrib.polyg(a1,a2,linepick).lines{destClust},'UserData',[a1 a2 linepick destClust]);
set(graphattrib.polyg(a1,a2,linepick).lines{destClust},'Color',figattrib.mixcolor(destClust,:));
set(graphattrib.polyg(a1,a2,linepick).lines{destClust},'UIContextMenu',figattrib.handles.polycontext);
graphattrib.polyg(a1,a2,linepick).type{destClust} = 1;

try
    clustattrib.clusters{destClust}.defineaxes = setdiff(clustattrib.clusters{clustnum}.defineaxes,[a1 a2 linepick],'rows');
    clustattrib.clusters{destClust}.defineaxes(end+1,:) = [a1 a2 linepick];
catch
    clustattrib.clusters{destClust}.defineaxes = [a1 a2 linepick];
end

findInPoints(a1,a2,linepick,destClust);
%-----------------------------------------------
function epochcopy(clustnum,source,dest,handles)
% allows the user to copy all polygons for cluster clustnum from time
% filter source to time filter dest.

global clustattrib;
global graphattrib;
global figattrib;
global clustdata;


def = clustattrib.clusters{clustnum}.defineaxes; 

releasepolydraw(handles);
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
        destfilt = [dest sourcefilt(2:end)];
        if (sourcefilt(1) == source)
                       
            linepick = -1;
       
            for j = 1:length(clustattrib.filterindex)
               if (isequal(clustattrib.filterindex{j},destfilt))
                    linepick = j;
                    break;
               end
            end
            if (linepick == -1)  %destination filter does not exist

                for j = 1:length(clustattrib.filterindex)
                    if isempty(clustattrib.filterindex{j})
                        linepick =  j;
                        break;
                    end
                end
            end
            if (linepick == -1)
                linepick = length(clustattrib.filterindex)+1;
            end

            clustattrib.filterindex{linepick} = destfilt;
            
            graphattrib.polyg(a1,a2,linepick).vertices{clustnum} = graphattrib.polyg(a1,a2,findex).vertices{clustnum};
            graphattrib.polyg(a1,a2,linepick).relativevertices{clustnum} = graphattrib.polyg(a1,a2,findex).relativevertices{clustnum};
            graphattrib.polyg(a1,a2,linepick).lines{clustnum} = line(graphattrib.polyg(a1,a2,findex).vertices{clustnum}(:,1),graphattrib.polyg(a1,a2,findex).vertices{clustnum}(:,2));
            set(graphattrib.polyg(a1,a2,linepick).lines{clustnum},'ButtonDownFcn','matclust(''polygon_ButtonDownFcn'',gcbo,guidata(gcbo))');

            set(graphattrib.polyg(a1,a2,linepick).lines{clustnum},'UserData',[a1 a2 linepick clustnum]);
            set(graphattrib.polyg(a1,a2,linepick).lines{clustnum},'Color',colorclip(figattrib.mixcolor(clustnum,:)-[.2 .2 .2]));
            set(graphattrib.polyg(a1,a2,linepick).lines{clustnum},'UIContextMenu',handles.polycontext);
            graphattrib.polyg(a1,a2,linepick).type{clustnum} = 1;
            try
                clustattrib.clusters{clustnum}.defineaxes = setdiff(clustattrib.clusters{clustnum}.defineaxes,[a1 a2 linepick],'rows');
                clustattrib.clusters{clustnum}.defineaxes(end+1,:) = [a1 a2 linepick];
            catch
                clustattrib.clusters{clustnum}.defineaxes = [a1 a2 linepick];
            end
            
            tmpfilter = fastandbit(clustdata.timefilters,destfilt(1));
            tmpfilter2 = fastandbit(clustdata.otherfilters,destfilt(2:end));
            clustdata.filteredpoints = tmpfilter & tmpfilter2;
            
            findInPoints(a1,a2,linepick,clustnum);
           
            curra1 = get(handles.listbox1,'UserData');
            curra2 = get(handles.listbox2,'UserData');
            [currfindex, filterline] = findfilterindex(handles);
            %hide the currently displayed polygons

            if ~isequal([curra1 curra2 currfindex],[a1 a2 linepick])
                set(graphattrib.polyg(a1,a2,linepick).lines{clustnum},'Visible','off');
            end          
            
        end
        
    end
    
    FilterPoints;
    updateclustinfo(handles);
    plotgraph(handles);
    addnewstate('copy polygons',handles); 
end
%------------------------------------------------
function turnoncluster(clust,handles,KeepCurrentHiddenStates)

global figattrib
global clustattrib

if (nargin < 3)
    KeepCurrentHiddenStates = 0;
end
set(figattrib.clustcontrol(clust,3:end),'Visible','on'); %turn on cluster control
set(figattrib.clustcontrol(clust,1),'BackgroundColor',figattrib.mixcolor(clust,:));

if (~KeepCurrentHiddenStates)
    set(figattrib.clustcontrol(clust,3),'Value',0); %turn filter off
    set(figattrib.clustcontrol(clust,4),'Value',1); %turn visible on
    clustattrib.hiddenclusters(clust,1:2) = 0;
else
    set(figattrib.clustcontrol(clust,3),'Value',clustattrib.hiddenclusters(clust,1)); 
    set(figattrib.clustcontrol(clust,4),'Value',1-clustattrib.hiddenclusters(clust,2)); 
end
if ~sum(ismember(clustattrib.clustersOn,clust))    
    clustattrib.clustersOn = [clustattrib.clustersOn;clust]; %add to 'on' list
end

eval(['set(handles.clustcopyMenu',num2str(clust),',''Enable'',''on'');'])
eval(['set(handles.epochcopyMenu',num2str(clust),',''Enable'',''on'');'])
eval(['set(handles.clustdeleteMenu',num2str(clust),',''Enable'',''on'');'])
eval(['set(handles.colorMenu',num2str(clust),',''Enable'',''on'');'])
eval(['set(handles.orderMenu',num2str(clust),',''Enable'',''on'');'])
eval(['set(handles.ClustToolsMenu',num2str(clust),',''Enable'',''on'');'])
%eval(['set(handles.viewISIMenu',num2str(clust),',''Enable'',''on'');'])
%eval(['set(handles.viewEvents',num2str(clust),',''Enable'',''on'');'])
%---------------------------------------------------
function turnoffcluster(clust,handles)

global figattrib
global clustattrib


set(figattrib.clustcontrol(clust,3:end),'Visible','off'); %turn off cluster control
set(figattrib.clustcontrol(clust,3),'Value',0); %turn filter off
set(figattrib.clustcontrol(clust,4),'Value',1); %turn visible on
set(figattrib.clustcontrol(clust,1),'BackgroundColor',figattrib.figureColor);
clustattrib.hiddenclusters(clust,1:2) = 0;
%remove all dependencies from the cluster
dep = find(clustattrib.dependencies(:,clust));
clustattrib.dependencies(clust,:) = 0;
clustattrib.dependencies(:,clust) = 0;
clustattrib.clustersOn = setdiff(clustattrib.clustersOn,clust,'rows'); %remove from 'on' list

eval(['set(handles.clustcopyMenu',num2str(clust),',''Enable'',''off'');'])
eval(['set(handles.epochcopyMenu',num2str(clust),',''Enable'',''off'');'])
eval(['set(handles.clustdeleteMenu',num2str(clust),',''Enable'',''off'');'])
eval(['set(handles.colorMenu',num2str(clust),',''Enable'',''off'');'])
eval(['set(handles.orderMenu',num2str(clust),',''Enable'',''off'');'])
eval(['set(handles.ClustToolsMenu',num2str(clust),',''Enable'',''off'');'])
%eval(['set(handles.viewISIMenu',num2str(clust),',''Enable'',''off'');'])
%eval(['set(handles.viewEvents',num2str(clust),',''Enable'',''off'');'])
set(figattrib.clustcontrol(clust,2),'TooltipString',['Cluster ',num2str(clust)]);
for i = 1:length(dep)
    findInPoints([],[],[],dep(i), find(clustattrib.dependencies(dep(i),:)));
end

%------------------------------------------------------
function paramrotate(hObject,handles)

global clustdata;
global graphattrib;
global figattrib;
global clustattrib;

try
    delete(figattrib.rotationwindowHandle);
end
h = rotationfunc(handles);
figattrib.rotationwindowHandle = h;

%---------------------------------------------------
function paramadd(hObject,handles)

%This function is currently not connected to the menu, as I found little
%use for adding new parameters.  The rotation function has taken its place.

global clustdata;
global graphattrib;
global figattrib;
global clustattrib;

origpath = pwd;
cd(figattrib.foldername);
cd ('ParameterPrograms');

if (ispc)
    [filename, pathname] = uigetfile('*.m','Open a parameter file');
else
    [filename, pathname] = filebrowse('open','filter','*.m','title','Open a parameter file');
end


eval(['cd ',pathname]);
try
    
    pointfind = strfind(filename,'.');
    filename = filename(1:pointfind-1);
    [outparam, newnames] = feval(filename);
catch
    if isempty(filename)
        return
    end
    error('Could not run m-file.  No paramters were added.')
    return        
end
cd(origpath);

if ~isequal(size(outparam,1),size(clustdata.params,1))
    error('Output of parameter generator is not of correct length')
    return    
end
if ~isequal(size(outparam,2),length(newnames))
    error('Output error in parametor genorator: the names output must be the same length as the number of colums of data')
    return
end

S = get(handles.listbox1,'String');
for i = 1:length(newnames)
    tmpfill = [clustdata.filledparam 0];
    firstzero = min(find(tmpfill == 0));
    clustdata.params(:,firstzero) = outparam(:,i);
    clustdata.names{firstzero} = newnames{i};
    clustdata.filledparam(firstzero) = 1;
    S{firstzero} = [num2str(firstzero),'  ',newnames{i}];
    clustdata.datarange(:,firstzero) = [min(clustdata.params(:,firstzero));max(clustdata.params(:,firstzero))]; %stores the minimum and maximum values (-/+ 10% of range) of each parameter
    clustdata.datarange(:,firstzero) = [min(clustdata.params(:,firstzero))-(.1*diff(clustdata.datarange(:,firstzero)));max(clustdata.params(:,firstzero))+(.1*diff(clustdata.datarange(:,firstzero)))];
    graphattrib.viewbox(firstzero,1:4,1:size(clustdata.params,2)) = repmat([0 1 0 1],[1 1 size(clustdata.params,2)]);
    graphattrib.viewbox(:,1:4,firstzero) = repmat([0 1 0 1],[size(clustdata.params,2) 1 1]);
end
set(handles.listbox1,'String',S);
set(handles.listbox2,'String',S);

%update the parameter file in the 'currently open' folder
%cd([figattrib.datafoldername,'M_dataopen']);
cd(figattrib.datafoldername);
cd('M_dataopen');
cd(figattrib.openfiles{figattrib.currentopenfile});
paramdata = clustdata.params;
save paramdata paramdata;
eval(['cd ',origpath]);

addnewstate('add parameter',handles); 

%------------------------------------------------------
function paramdelete(hObject,handles)

%This function is currently not connected to the menu, as I found little
%use for adding new parameters.  The rotation function has taken its place.

global clustdata;
global graphattrib;
global clustattrib;

listnum = get(hObject,'UserData');

eval(['currval = get(handles.listbox',num2str(listnum),',''Value'');'])

if (currval>clustdata.origparams)
    reply = questdlg('Warning: this will delete all cluster boxes using this parameter. It cannot be undone!','Warning','Cancel','Delete parameter','Delete parameter');
    if strcmp(reply, 'Delete parameter')
        releasepolydraw(handles);
        S = get(handles.listbox1,'String');
        S{currval} = '';
        clustdata.filledparam(currval) = 0;
        clustdata.names{currval} = '';

        
        %memoryindex = find(clustdata.timefiltermemmap == currval);
	
	    %newval = max(clustdata.timefiltermemmap(find(clustdata.timefiltermemmap<currval)));
	
		a1 = get(handles.listbox1,'UserData');
		a2 = get(handles.listbox2,'UserData');
		[findex, filterline] = findfilterindex(handles);
		%hide the currently displayed polygons
		try
            for i = 1:length(graphattrib.polyg(a1,a2,findex).lines)
                try
                    set(graphattrib.polyg(a1,a2,findex).lines{i},'Visible','off');
                    delete(graphattrib.polyg(a1,a2,findex).highlight{i});    
                end
            end
		end
	
	    set(handles.listbox1,'String',S);
        set(handles.listbox2,'String',S);
	    set(handles.listbox1,'Value',1);
        set(handles.listbox2,'Value',2);
        set(handles.listbox1,'UserData',1);
        set(handles.listbox2,'UserData',2);
        set(handles.paramdeleteContext1,'Enable','off');
        set(handles.paramdeleteContext2,'Enable','off');
        a1 = get(handles.listbox1,'UserData');
		a2 = get(handles.listbox2,'UserData');
	
		%show the new polygons in in the now active filter
		try
            for i = 1:length(graphattrib.polyg(a1,a2,findex).lines)
                try
                    set(graphattrib.polyg(a1,a2,findex).lines{i},'Visible','on');        
                end
            end
		end    
	
	    for c = clustattrib.clustersOn'
            dpoly = min(find((clustattrib.clusters{c}.defineaxes(:,1) == currval)|(clustattrib.clusters{c}.defineaxes(:,2) == currval)));
            while ~isempty(dpoly)
                deletepoly(clustattrib.clusters{c}.defineaxes(dpoly,1),clustattrib.clusters{c}.defineaxes(dpoly,2),clustattrib.clusters{c}.defineaxes(dpoly,3),c);
                dpoly = min(find((clustattrib.clusters{c}.defineaxes(:,1) == currval)|(clustattrib.clusters{c}.defineaxes(:,2) == currval)));    
            end
        end
	
	
	    updateclustinfo(handles);
	
	    plotgraph(handles);
        addnewstate('delete parameter',handles); 
    end
end
%---------------------------------------------------
function paramrotationdelete(hObject,handles)

global clustdata;
global graphattrib;
global clustattrib;

listnum = get(hObject,'UserData');

eval(['currval = get(handles.listbox',num2str(listnum),',''Value'');'])

if (currval>clustdata.origparams)
    reply = questdlg('Warning: this will delete all cluster boxes using this rotation.','Warning','Cancel','Delete rotation','Delete rotation');
    if strcmp(reply, 'Delete rotation')
        releasepolydraw(handles);
        S = get(handles.listbox1,'String');
        S{currval} = '';
        clustdata.filledparam(currval) = 0;
        clustdata.names{currval} = '';

        
        %memoryindex = find(clustdata.timefiltermemmap == currval);
	
	    %newval = max(clustdata.timefiltermemmap(find(clustdata.timefiltermemmap<currval)));
	
		a1 = get(handles.listbox1,'UserData');
		a2 = get(handles.listbox2,'UserData');
		[findex, filterline] = findfilterindex(handles);
		%hide the currently displayed polygons
		try
            for i = 1:length(graphattrib.polyg(a1,a2,findex).lines)
                try
                    set(graphattrib.polyg(a1,a2,findex).lines{i},'Visible','off');
                    delete(graphattrib.polyg(a1,a2,findex).highlight{i});    
                end
            end
		end
	
	    set(handles.listbox1,'String',S);
        set(handles.listbox2,'String',S);
	    set(handles.listbox1,'Value',1);
        set(handles.listbox2,'Value',2);
        set(handles.listbox1,'UserData',1);
        set(handles.listbox2,'UserData',2);
        set(handles.paramdeleteContext1,'Enable','off');
        set(handles.paramdeleteContext2,'Enable','off');
        a1 = get(handles.listbox1,'UserData');
		a2 = get(handles.listbox2,'UserData');
	
		%show the new polygons in in the now active filter
		try
            for i = 1:length(graphattrib.polyg(a1,a2,findex).lines)
                try
                    set(graphattrib.polyg(a1,a2,findex).lines{i},'Visible','on');        
                end
            end
		end    
	
	    for c = clustattrib.clustersOn'
            dpoly = min(find((clustattrib.clusters{c}.defineaxes(:,1) == currval)|(clustattrib.clusters{c}.defineaxes(:,2) == currval)));
            while ~isempty(dpoly)
                deletepoly(clustattrib.clusters{c}.defineaxes(dpoly,1),clustattrib.clusters{c}.defineaxes(dpoly,2),clustattrib.clusters{c}.defineaxes(dpoly,3),c);
                dpoly = min(find((clustattrib.clusters{c}.defineaxes(:,1) == currval)|(clustattrib.clusters{c}.defineaxes(:,2) == currval)));    
            end
        end
	
	
	    updateclustinfo(handles);
	
	    plotgraph(handles);
        addnewstate('delete rotation',handles); 
    end
end
%------------------------------------------------
function parameditname(hObject,handles)

global clustdata;

a = inputdlg('Enter new parameter name','Edit name');

if ~isempty(a)
    listnum = get(hObject,'UserData');
    eval(['currval = get(handles.listbox',num2str(listnum),',''Value'');'])
    
    clustdata.names{currval} = [num2str(currval),'  ',a{1}];
    set(handles.listbox1,'String',clustdata.names);
    set(handles.listbox2,'String',clustdata.names);
    addnewstate('edit name',handles); 
end
%-----------------------------------------------------
function zoomFull()

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;
handles = figattrib.handles;

releasepolydraw(handles);
currx = get(handles.listbox1,'Value');
curry = get(handles.listbox2,'Value');
graphattrib.viewbox(currx,1:4,curry) = [0 1 0 1];   
plotgraph(handles)
addnewstate('zoom',handles); 
%------------------------------------------------
function zoomIn(cursorLoc)

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;
handles = figattrib.handles;


releasepolydraw(handles);
currx = get(handles.listbox1,'Value');
curry = get(handles.listbox2,'Value');
currview = graphattrib.viewbox(currx,1:4,curry);
xrange = currview(2)-currview(1);
yrange = currview(4)-currview(3);
if (nargin > 0)
    viewrange = get(graphattrib.graphwindow,'UserData');
    relativePoint = [(cursorLoc(1,1)+viewrange(3))/viewrange(1) (cursorLoc(1,2)+viewrange(5))/viewrange(2)];
    diffToCursor = [relativePoint(1) - (currview(1)+(xrange/2)) relativePoint(2)-(currview(3)+(yrange/2))];
else
    diffToCursor = [0 0];
end
%currview = [(currview(1)+.05*(xrange)) (currview(2)-.05*(xrange)) (currview(3)+.05*(yrange)) (currview(4)-.05*(yrange))];
panMult = .1;
currview = [(currview(1)+.05*(xrange))+panMult*diffToCursor(1) (currview(2)-.05*(xrange))+panMult*diffToCursor(1) (currview(3)+.05*(yrange))+panMult*diffToCursor(2) (currview(4)-.05*(yrange))+panMult*diffToCursor(2)];

%currview = [min([currview(2)-.1 currview(1)+.1]) max([currview(1)+.1 currview(2)-.1]) min([currview(4)-.1 currview(3)+.1]) max([currview(3)+.1 currview(4)-.1])];
graphattrib.viewbox(currx,1:4,curry) = currview;   
plotgraph(handles)
addnewstate('zoom',handles); 
%---------------------------------------------------
function zoomAllIn()

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;
handles = figattrib.handles;

releasepolydraw(handles);
for currx = 1:size(clustdata.params,2)
    for curry = 1:size(clustdata.params,2)
        
        currview = graphattrib.viewbox(currx,1:4,curry);
        xrange = currview(2)-currview(1);
        yrange = currview(4)-currview(3);
        currview = [(currview(1)+.25*(xrange)) (currview(2)-.25*(xrange)) (currview(3)+.25*(yrange)) (currview(4)-.25*(yrange))];
        graphattrib.viewbox(currx,1:4,curry) = currview;   
    end
end
plotgraph(handles)
addnewstate('zoom',handles); 
%------------------------------------------------
function zoomOut(cursorLoc)

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;
handles = figattrib.handles;

releasepolydraw(handles);
currx = get(handles.listbox1,'Value');
curry = get(handles.listbox2,'Value');
currview = graphattrib.viewbox(currx,1:4,curry);
xrange = currview(2)-currview(1);
yrange = currview(4)-currview(3);

if (nargin > 0)
    viewrange = get(graphattrib.graphwindow,'UserData');
    relativePoint = [(cursorLoc(1,1)+viewrange(3))/viewrange(1) (cursorLoc(1,2)+viewrange(5))/viewrange(2)];
    diffToCursor = [relativePoint(1) - (currview(1)+(xrange/2)) relativePoint(2)-(currview(3)+(yrange/2))];
else
    diffToCursor = [0 0];
end
%currview = [(currview(1)+.05*(xrange)) (currview(2)-.05*(xrange)) (currview(3)+.05*(yrange)) (currview(4)-.05*(yrange))];
panMult = .1;
currview = [max([0 currview(1)-.05-panMult*diffToCursor(1)]) min([1 currview(2)+.05-panMult*diffToCursor(1)]) max([0 currview(3)-.05-panMult*diffToCursor(2)]) min([1 currview(4)+.05-panMult*diffToCursor(2)])];
graphattrib.viewbox(currx,1:4,curry) = currview;   
plotgraph(handles)
addnewstate('zoom',handles); 
%------------------------------------------------
function zoomAllOut()

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;
handles = figattrib.handles;
releasepolydraw(handles);

for currx = 1:size(clustdata.params,2)
    for curry = 1:size(clustdata.params,2)
        currview = graphattrib.viewbox(currx,1:4,curry);
        currview = [max([0 currview(1)-.05]) min([1 currview(2)+.05]) max([0 currview(3)-.05]) min([1 currview(4)+.05])];
        graphattrib.viewbox(currx,1:4,curry) = currview;
    end
end
plotgraph(handles)
addnewstate('zoom',handles); 
%--------------------------------------------------
function setZoomRange(percentile)
%function to zoom to a particular data percentage (100-percentage and
%percentage of data)

global clustdata;
global clustattrib
global graphattrib;
global figattrib;
handles = figattrib.handles;

prcshow = percentile;
%clustdata.datarange = [min(clustdata.params)-(.1*diff(clustdata.datarange));max(clustdata.params)+(.1*diff(clustdata.datarange))];

for i = 1:size(clustdata.params,2)
    percentages = prctile(clustdata.params(:,i),[0 100-prcshow prcshow 100]);
    datarange = clustdata.datarange(:,i);
    viewboxvalues1 = (percentages(2)-datarange(1))/(datarange(2)-datarange(1));
    viewboxvalues2 = (percentages(3)-datarange(1))/(datarange(2)-datarange(1));
   
    graphattrib.viewbox(i,1,:) = viewboxvalues1;
    graphattrib.viewbox(i,2,:) = viewboxvalues2;
    graphattrib.viewbox(:,3,i) = viewboxvalues1;
    graphattrib.viewbox(:,4,i) = viewboxvalues2;   
end

plotgraph(handles);
%---------------------------------------------------------------------
function setPanelWidth(width,saveVal)
if (nargin < 2)
    saveVal = 1; %the default is to save the new panel width in the user defaults file
end

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;
handles = figattrib.handles;

if (width < 265)
    width = 265;
end
figattrib.sidePanelWidth = width; %the width of the side control panel
if (saveVal)
    currdir = pwd;
    cd(figattrib.foldername);
    load matclust_defaults;
    matclust_defaults.sidePanelWidth = width;
    save('matclust_defaults','matclust_defaults');
    cd(currdir);
end
figure1_ResizeFcn(handles.figure1,handles)
    
%------------------------------------------------------
function setFilterBoxDivider(dividerLoc)

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;
handles = figattrib.handles;

if ((dividerLoc < .1)|(dividerLoc > .9))
    error('Filter box divider location must be between .1 and .9');
end
figattrib.filterBoxDivider = dividerLoc; %the division point between the two filter list boxes (range between 0 and 1) 
currdir = pwd;
cd(figattrib.foldername);
load matclust_defaults;
matclust_defaults.filterBoxDivider = dividerLoc;
save('matclust_defaults','matclust_defaults');
cd(currdir);
figure1_ResizeFcn(handles.figure1,handles)

%--------------------------------------------------
function releasepolydraw(handles)

global clustattrib;
global graphattrib;
global clustdata;

if (graphattrib.polydraw == 1)
	a1 = get(handles.listbox1,'UserData');
	a2 = get(handles.listbox2,'UserData');
	[findex, filterline] = findfilterindex(handles);
    findex = findex(1);
	graphattrib.polydraw = 0;
	findInPoints(a1,a2,findex,clustattrib.currclust);  
    if (length(find(clustdata.timefiltersOn)) > 1)
        tfilts = find(clustdata.timefiltersOn);
        for tnum = 2:length(tfilts)
            copypoly(a1,a2,findex,clustattrib.currclust,tfilts(tnum)); %copy the polygon
        end
    end
    updateclustinfo(handles);
	plotgraph(handles);  
    highlightCurrentCluster(handles);
    addnewstate('draw polygon',handles); 
    
end
%-----------------------------------------------------
function openrawparam_Callback(hObject, handles)
global clustdata;
global graphattrib;
global clustattrib;
global figattrib;    

tmpclustdata = clustdata;
tmpgraphattrib = graphattrib;
tmpclustattrib = clustattrib;
tmpfigattrib = figattrib;

if (ispc)
    [filename, pathname] = uigetfile('*.mat','Open a parameter file');
else
    [filename, pathname] = filebrowse('open','filter','*.mat','title','Open a parameter file');
end

if ((filename))
    if ~(clustattrib.nodata)
        savefigstate
    end
    origpath = pwd;
	cd(pathname);
    set(figattrib.handles.statusText,'String','Loading data...');
    set(figattrib.handles.statusText,'UserData','Loading data...');
    drawnow;
	
	load(filename);
    %graphattrib.graphwindow = tmpgraphattrib.graphwindow; 
    %graphattrib.Himage = tmpgraphattrib.Himage;
    %cd(origpath);
   
    try
        clearworkspace(1,1);
		clearhistory;
        if ((exist('filedata'))&(isstruct(filedata)))
            if isfield(filedata,'filename')
                clustattrib.datafile = filedata.filename;
            end
             if isfield(filedata,'params')
                clustdata.params = filedata.params;            
                  if isfield(filedata,'paramnames')
                    clustdata.names = filedata.paramnames;
                  else
                      for i = 1:size(clustdata.params,2)
                          clustdata.names{i,1} = [num2str(i), ' Column ',num2str(i)];
                      end
                  end
             else
                error('Input is not valid');
             end
             if isfield(filedata,'customvar')
                 clustdata.customvar = filedata.customvar;
             end
        elseif ((exist('filedata'))&(isnumeric(filedata)))
            clustdata.params = filedata;
            for i = 1:size(clustdata.params,2)
                 clustdata.names{i,1} = [num2str(i), ' Column ',num2str(i)];
            end
         else
            error('File is not valid');   
        end
        %clustattrib.datafileloc = filedata.fileloc;
		clustdata.origparams = length(clustdata.names);
		for i = 1:length(clustdata.names)
            clustdata.names{i} = [num2str(i),'  ',clustdata.names{i}];
            clustdata.filledparam(i) = 1;
		end
		
		clustdata.timefilters = int32(zeros(size(clustdata.params,1),1)); %contains the filtere points of up to 32 time filters
		clustdata.timefilters = fastbitset(clustdata.timefilters,1,logical(1)); %the first filter is all times
		clustdata.timefiltermemmap = zeros(32,1); %because the flter order can be changed by the user, this map translates the filter numbers to those stored in memory 
		clustdata.timefiltermemmap(1) = 1;
		clustdata.timefiltersOn = clustdata.timefiltermemmap'; %which time filter is currently on
		clustdata.otherfilters = int32(zeros(size(clustdata.params,1),1));
		clustdata.otherfilters = fastbitset(clustdata.otherfilters,1,logical(1));
		clustdata.otherfiltermemmap = zeros(32,1);
		clustdata.otherfiltermemmap(1) = 1;
		clustdata.otherfiltersOn = clustdata.otherfiltermemmap'; %for the 'other' filters, more than one filter can be on at the same time 
		clustdata.filtermemmap = [clustdata.timefiltermemmap; clustdata.otherfiltermemmap];
		clustdata.filteredpoints = logical(ones(size(clustdata.params,1),1)); %stores which points are currently let through the filters
		clustdata.datarange = [min(clustdata.params);max(clustdata.params)]; %stores the minimum and maximum values (-/+ 10% of range) of each parameter
		clustdata.datarange = [min(clustdata.params)-(.1*diff(clustdata.datarange));max(clustdata.params)+(.1*diff(clustdata.datarange))];
		baddim = find(diff(clustdata.datarange)==0);
        clustdata.datarange(1:2,baddim) = [clustdata.datarange(1,baddim)-1;clustdata.datarange(2,baddim)+1]; 
        graphattrib.viewbox = repmat([0 1 0 1],[size(clustdata.params,2) 1 size(clustdata.params,2)]);  %stores the current view window for each axis pair [xmin xmax ymin ymax]
		graphattrib.oldviewbox = graphattrib.viewbox;
		clustattrib.pointexclude = int32(zeros(size(clustdata.params,1),1)); %contains double-precision vectors storing the outside-inside status of each data point for 32 different polygons
		clustattrib.pointinclude = int32(zeros(size(clustdata.params,1),1));
		clustattrib.dependencies = false(size(figattrib.mixcolor,1));
        
		%clustattrib.cluster0 = [1:size(clustdata.params,1)]';
		set(handles.listbox1,'Value',1);
		set(handles.listbox1,'UserData',1);
		set(handles.listbox1,'String',clustdata.names);
		set(handles.listbox2,'Value',2);
		set(handles.listbox2,'UserData',2);
		set(handles.listbox2,'String',clustdata.names);
		set(handles.paramdeleteContext1,'Enable','off');
		set(handles.paramdeleteContext2,'Enable','off');
		clustdata.timefilterranges = clustdata.datarange(:,1)';
		clustattrib.cluster0attrib.show = 1;
        set(handles.cluster0button,'String','Cluster 0: visible');	
        set(handles.cluster0button,'Value',1);	
	
        clustattrib.nodata = 0;
        set(handles.file_closeMenu,'Enable','on');
        set(handles.file_saveasMenu,'Enable','on');
        set(handles.file_saveMenu,'Enable','on');
        set(handles.file_curropenMenu,'Enable','on');
        set(handles.file_exportMenu,'Enable','on');
        set(handles.listbox1,'Enable','on');
        set(handles.listbox2,'Enable','on');
        set(handles.timefilterList,'Enable','on');
        set(handles.otherfilterList,'Enable','on');
        set(handles.paramaddContext1,'Enable','on');
        set(handles.paramaddContext2,'Enable','on');
        set(handles.parameditContext1,'Enable','on');
        set(handles.parameditContext2,'Enable','on');
        set(handles.paramdeleteContext1,'Enable','off');
        set(handles.paramdeleteContext2,'Enable','off');
        clustattrib.currentparamfilename = filename;
        clustattrib.currentfilepath = pathname;
		clustattrib.newchanges = 1;
        open = OpenNewDataFolder;
        
        currdir = pwd;
        %cd([figattrib.datafoldername,'M_dataopen']);
        cd(figattrib.datafoldername);
        cd('M_dataopen');
        cd(figattrib.openfiles{figattrib.currentopenfile});
        currstate = 1;
        clustattrib.currstate = currstate;
        %maxstate = 1;
        paramdata = clustdata.params;
        save paramdata paramdata;
        %clustdata = rmfield(clustdata,'params');
        clustattrib.states{1}.clustdata = rmfield(clustdata,'params');
        clustattrib.states{1}.clustattrib = rmfield(clustattrib,'states');
        clustattrib.states{1}.graphattrib = graphattrib;
        states = clustattrib.states;
        save states states;
        %save state1 clustdata graphattrib clustattrib;
        %clustdata.params = paramdata;
        save currstate currstate;
        %save maxstate maxstate;
        cd(currdir);
        set(handles.undoMenu,'Label','Undo');
        set(handles.undoMenu,'Enable','off'); 
        set(handles.redoMenu,'Enable','off'); 
        set(handles.figure1,'Name',['MATCLUST  ',figattrib.openfiles{figattrib.currentopenfile}]);
        updatefileswitcher(handles);
        
        plotgraph(handles); 
    catch  
        clustdata = tmpclustdata;
        graphattrib = tmpgraphattrib;
        clustattrib = tmpclustattrib;
        figattrib = tmpfigattrib;
        set(figattrib.handles.statusText,'String','Ready');
        set(figattrib.handles.statusText,'UserData','Ready');

        error('Not a correct parameter file')
   end
   set(figattrib.handles.statusText,'String','Ready');
   set(figattrib.handles.statusText,'UserData','Ready');
    
end
%----------------------------------------------------
function openinputparam_Callback(handles,filedata)
global clustdata;
global graphattrib;
global clustattrib;
global figattrib;    

tmpclustdata = clustdata;
tmpgraphattrib = graphattrib;
tmpclustattrib = clustattrib;
tmpfigattrib = figattrib;



if (1)
    if ~(clustattrib.nodata)
        savefigstate
    end
    origpath = pwd;
	
    %graphattrib.graphwindow = tmpgraphattrib.graphwindow; 
    %graphattrib.Himage = tmpgraphattrib.Himage;
    %cd(origpath);
   
    try
        clearworkspace(1,1);
		clearhistory;
        if ((exist('filedata'))&(isstruct(filedata)))
            if isfield(filedata,'filename')
                clustattrib.datafile = filedata.filename;
            end
             if isfield(filedata,'params')
                clustdata.params = filedata.params;            
                  if isfield(filedata,'paramnames')
                    clustdata.names = filedata.paramnames;
                  else
                      for i = 1:size(clustdata.params,2)
                          clustdata.names{i,1} = [num2str(i), ' Column ',num2str(i)];
                      end
                  end
             else
                error('Input is not valid');
             end
             if isfield(filedata,'customvar')
                 clustdata.customvar = filedata.customvar;
             end
        elseif ((exist('filedata'))&(isnumeric(filedata)))
            clustdata.params = filedata;
            for i = 1:size(clustdata.params,2)
                 clustdata.names{i,1} = [num2str(i), ' Column ',num2str(i)];
            end
        else
            error('Input is not valid');
        end
        %clustattrib.datafileloc = filedata.fileloc;
		clustdata.origparams = length(clustdata.names);
		for i = 1:length(clustdata.names)
            clustdata.names{i} = [num2str(i),'  ',clustdata.names{i}];
            clustdata.filledparam(i) = 1;
		end
		
		clustdata.timefilters = int32(zeros(size(clustdata.params,1),1)); %contains the filtere points of up to 32 time filters
		clustdata.timefilters = fastbitset(clustdata.timefilters,1,logical(1)); %the first filter is all times
		clustdata.timefiltermemmap = zeros(32,1); %because the flter order can be changed by the user, this map translates the filter numbers to those stored in memory 
		clustdata.timefiltermemmap(1) = 1;
		clustdata.timefiltersOn = clustdata.timefiltermemmap'; %which time filter is currently on
		clustdata.otherfilters = int32(zeros(size(clustdata.params,1),1));
		clustdata.otherfilters = fastbitset(clustdata.otherfilters,1,logical(1));
		clustdata.otherfiltermemmap = zeros(32,1);
		clustdata.otherfiltermemmap(1) = 1;
		clustdata.otherfiltersOn = clustdata.otherfiltermemmap'; %for the 'other' filters, more than one filter can be on at the same time 
		clustdata.filtermemmap = [clustdata.timefiltermemmap; clustdata.otherfiltermemmap];
		clustdata.filteredpoints = logical(ones(size(clustdata.params,1),1)); %stores which points are currently let through the filters
		clustdata.datarange = [min(clustdata.params);max(clustdata.params)]; %stores the minimum and maximum values (-/+ 10% of range) of each parameter
		clustdata.datarange = [min(clustdata.params)-(.1*diff(clustdata.datarange));max(clustdata.params)+(.1*diff(clustdata.datarange))];
		baddim = find(diff(clustdata.datarange)==0);
        clustdata.datarange(1:2,baddim) = [clustdata.datarange(1,baddim)-1;clustdata.datarange(2,baddim)+1]; 
        graphattrib.viewbox = repmat([0 1 0 1],[size(clustdata.params,2) 1 size(clustdata.params,2)]);  %stores the current view window for each axis pair [xmin xmax ymin ymax]
		graphattrib.oldviewbox = graphattrib.viewbox;
		clustattrib.pointexclude = int32(zeros(size(clustdata.params,1),1)); %contains double-precision vectors storing the outside-inside status of each data point for 32 different polygons
		clustattrib.pointinclude = int32(zeros(size(clustdata.params,1),1));
        clustattrib.dependencies = false(size(figattrib.mixcolor,1));
		
		%clustattrib.cluster0 = [1:size(clustdata.params,1)]';
		set(handles.listbox1,'Value',1);
		set(handles.listbox1,'UserData',1);
		set(handles.listbox1,'String',clustdata.names);
		set(handles.listbox2,'Value',2);
		set(handles.listbox2,'UserData',2);
		set(handles.listbox2,'String',clustdata.names);
		set(handles.paramdeleteContext1,'Enable','off');
		set(handles.paramdeleteContext2,'Enable','off');
		clustdata.timefilterranges = clustdata.datarange(:,1)';
		clustattrib.cluster0attrib.show = 1;
        set(handles.cluster0button,'String','Cluster 0: visible');	
        set(handles.cluster0button,'Value',1);	
	
        clustattrib.nodata = 0;
        set(handles.file_closeMenu,'Enable','on');
        set(handles.file_saveasMenu,'Enable','on');
        set(handles.file_saveMenu,'Enable','on');
        set(handles.file_curropenMenu,'Enable','on');
        set(handles.file_exportMenu,'Enable','on');
        set(handles.listbox1,'Enable','on');
        set(handles.listbox2,'Enable','on');
        set(handles.timefilterList,'Enable','on');
        set(handles.otherfilterList,'Enable','on');
        set(handles.paramaddContext1,'Enable','on');
        set(handles.paramaddContext2,'Enable','on');
        set(handles.parameditContext1,'Enable','on');
        set(handles.parameditContext2,'Enable','on');
        set(handles.paramdeleteContext1,'Enable','off');
        set(handles.paramdeleteContext2,'Enable','off');
        clustattrib.currentparamfilename = 'inputdata';
		clustattrib.newchanges = 1;
        open = OpenNewDataFolder;
        
        currdir = pwd;
        %cd([figattrib.datafoldername,'M_dataopen']);
        cd(figattrib.datafoldername);
        cd('M_dataopen');
        cd(figattrib.openfiles{figattrib.currentopenfile});
        currstate = 1;
        clustattrib.currstate = currstate;
        %maxstate = 1;
        paramdata = clustdata.params;
        save paramdata paramdata;
        %clustdata = rmfield(clustdata,'params');
        clustattrib.states{1}.clustdata = rmfield(clustdata,'params');
        clustattrib.states{1}.clustattrib = rmfield(clustattrib,'states');
        clustattrib.states{1}.graphattrib = graphattrib;
        states = clustattrib.states;
        save states states;
        %save state1 clustdata graphattrib clustattrib;
        %clustdata.params = paramdata;
        save currstate currstate;
        %save maxstate maxstate;
        cd(currdir);
        set(handles.undoMenu,'Label','Undo');
        set(handles.undoMenu,'Enable','off'); 
        set(handles.redoMenu,'Enable','off'); 
        set(handles.figure1,'Name',['MATCLUST  ',figattrib.openfiles{figattrib.currentopenfile}]);
        updatefileswitcher(handles);
        
        plotgraph(handles); 
    catch  
        clustdata = tmpclustdata;
        graphattrib = tmpgraphattrib;
        clustattrib = tmpclustattrib;
        figattrib = tmpfigattrib;
        disp('Error with input-- not a valid format');
   end
   set(figattrib.handles.statusText,'String','Ready');
   set(figattrib.handles.statusText,'UserData','Ready');
    
end
%----------------------------------------------
function closeMenu_Callback(hObject,handles)

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;

go = 1;
if (clustattrib.newchanges)
    if ~isempty(clustattrib.currentfilename)
        reply = questdlg(['Save changes to ',clustattrib.currentfilename,'?'],'Save changes?','Cancel','No','Save changes','Save changes');
    else
        reply = questdlg(['Save changes to ',figattrib.openfiles{figattrib.currentopenfile},'?'],'Save changes?','Cancel','No','Save changes','Save changes');
    end
    if (strcmp(reply,'Save changes'))
        saveMenu_Callback(handles, handles);
        if (clustattrib.newchanges == 0)
            go = 1;
        else
            go = 0;
        end
    elseif (strcmp(reply,'No'))
        go = 1;
    elseif (strcmp(reply,'Cancel'))
        go = 0;
    end
end
if (go)        
	currdir = pwd;
	%cd([figattrib.datafoldername,'M_dataopen']);
	cd(figattrib.datafoldername);
    cd('M_dataopen');
    if (ispc)
         system(['rd /S /Q ',figattrib.openfiles{figattrib.currentopenfile}]);
        
	elseif (isunix)
        system(['rm -r -f ', figattrib.openfiles{figattrib.currentopenfile}]);
	end
	figattrib.currentopenfile = [];
	dirnames = dir;
	figattrib.openfiles = [];
	count = 1;
	if (length(dirnames)>2)
        for i = 3:length(dirnames)
            if (dirnames(i).isdir)
                tmpname = dirnames(i).name;
                figattrib.openfiles{count} = tmpname;    
                count = count+1;    
            end
        end    
	end
	cd(currdir);
	if ~isempty(figattrib.openfiles)    
        figattrib.currentopenfile = 1;
        changefilefocus_after_close(figattrib.openfiles{figattrib.currentopenfile},1,1,1,0,handles)  
	else
        clearworkspace(1,1);
        clearhistory;
        set(handles.file_closeMenu,'Enable','off');
        set(handles.file_saveasMenu,'Enable','off');
        set(handles.file_saveMenu,'Enable','off');
        set(handles.file_curropenMenu,'Enable','off');
        set(handles.file_exportMenu,'Enable','off');
        set(handles.listbox1,'Enable','off');
        set(handles.listbox2,'Enable','off');
        set(handles.timefilterList,'Enable','off');
        set(handles.otherfilterList,'Enable','off');
        set(handles.paramaddContext1,'Enable','off');
        set(handles.paramaddContext2,'Enable','off');
        set(handles.parameditContext1,'Enable','off');
        set(handles.parameditContext2,'Enable','off');
        set(handles.paramdeleteContext1,'Enable','off');
        set(handles.paramdeleteContext2,'Enable','off');
        set(handles.redoMenu,'Enable','off');
        set(handles.undoMenu,'Label','Undo');
        set(handles.undoMenu,'Enable','off');
        setclustsizeinfo;
	end
	updatefileswitcher(handles);
    
end
%---------------------------------------------
function clearhistory()
global clustattrib;

clustattrib.states = [];
clustattrib.currstate = 1;

%------------------------------------------------
function clearworkspace(resetaxes,clearfilters)

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;

for clustnum = clustattrib.clustersOn' 
	def = clustattrib.clusters{clustnum}.defineaxes;	
	if ~isempty(def)       
        for i = 1:size(def,1)
            a1 = def(i,1);
            a2 = def(i,2);
            findex = def(i,3);
            try 
                delete(graphattrib.polyg(a1,a2,findex).lines{clustnum});
			end
			%turn highlighting off
			try 
                delete(graphattrib.polyg(a1,a2,findex).highlight{clustnum});             
			end
			%remove the record of a box in these axes if one exists
			
			graphattrib.polyg(a1,a2,findex).vertices{clustnum} = [];
			graphattrib.polyg(a1,a2,findex).relativevertices{clustnum} = [];
			clustattrib.clusters{clustnum}.box{a1,a2,findex} = [];

        end
        
    end
    clustattrib.clusters{clustnum}.defineaxes = [];    
    clustattrib.clusters{clustnum}.index = [];
    clustattrib.clusters{clustnum}.polyindex = [];
    turnoffcluster(clustnum,figattrib.handles);
end
    
clustattrib.clusters = [];
clustattrib.filterindex = [];
graphattrib.polyg = [];
graphattrib.relativepolypress = [];
%for i = clustattrib.clustersOn' 
%    clustdelete(i, figattrib.handles, 0);
    
%end

clustattrib.currentfilepath = [];
clustattrib.currentfilename = []; %stores the name of the current MatClust file
clustattrib.currentparamfilename = [];

handles = figattrib.handles;
clustdata.params = [1 1;2 2];
clustdata.origparams = 1;
clustdata.names = {'',
                    ' '};

clustdata.filledparam = [];
clustdata.timefilterranges = [];
clustattrib.lastaction = [];
clustdata.datarange = [min(clustdata.params,[],1);max(clustdata.params,[],1)]; %stores the minimum and maximum values (-/+ 10% of range) of each parameter
%clustdata.datarange = [min(clustdata.params)-(.1*diff(clustdata.datarange));max(clustdata.params)+(.1*diff(clustdata.datarange))];
graphattrib.viewbox = repmat([0 1 0 1],[size(clustdata.params,2) 1 size(clustdata.params,2)]);  %stores the current view window for each axis pair [xmin xmax ymin ymax]
graphattrib.oldviewbox = graphattrib.viewbox;
clustattrib.takenpolys = [];
clustattrib.pointexclude = int32(zeros(size(clustdata.params,1),1)); %contains 32bit-precision vectors storing the outside-inside status of each data point for 32 different polygons
clustattrib.pointinclude = int32(zeros(size(clustdata.params,1),1));
clustattrib.eventeditindex = [];
%clustattrib.cluster0 = [1:size(clustdata.params,1)]';

if (resetaxes)
	set(handles.listbox1,'Value',1);
	set(handles.listbox1,'UserData',1);
	set(handles.listbox1,'String',clustdata.names);
	set(handles.listbox2,'Value',1);
	set(handles.listbox2,'UserData',1);
	set(handles.listbox2,'String',clustdata.names);
	set(handles.paramdeleteContext1,'Enable','off');
	set(handles.paramdeleteContext2,'Enable','off');
end   
 if (clearfilters)   
	for bfill = 1:32
        blanks(bfill) = {[num2str(bfill)]};
	end
	blanks{1} = [num2str(1),187,' All points'];
	clustdata.timefilternames = blanks;
	clustdata.otherfilternames = blanks;
	clustdata.timefilterranges = clustdata.datarange(:,1)';
	set(handles.timefilterList,'String',blanks);
	set(handles.otherfilterList,'String',blanks);
	set(handles.timefilterList,'UserData',[1 0]);
	set(handles.otherfilterList,'UserData',[1 0]);
    clustdata.timefilters = int32(zeros(size(clustdata.params,1),1)); %contains the filtere points of up to 32 time filters
	clustdata.timefilters = fastbitset(clustdata.timefilters,1,logical(1)); %the first filter is all times
	clustdata.timefiltermemmap = zeros(32,1); %because the flter order can be changed by the user, this map translates the filter numbers to those stored in memory 
	clustdata.timefiltermemmap(1) = 1;
	clustdata.timefiltersOn = clustdata.timefiltermemmap'; %which time filter is currently on
	clustdata.otherfilters = int32(zeros(size(clustdata.params,1),1));
	clustdata.otherfilters = fastbitset(clustdata.otherfilters,1,logical(1));
	clustdata.otherfiltermemmap = zeros(32,1);
	clustdata.otherfiltermemmap(1) = 1;
	clustdata.otherfiltersOn = clustdata.otherfiltermemmap'; %for the 'other' filters, more than one filter can be on at the same time 
	clustdata.filtermemmap = [clustdata.timefiltermemmap; clustdata.otherfiltermemmap];
	clustdata.filteredpoints = logical(ones(size(clustdata.params,1),1)); %stores which points are currently let through the filters

    resetfilters(0);
end

updateclustinfo(figattrib.handles);
set(handles.figure1,'Name','MATCLUST');

set(handles.mainClustMenu,'Enable','off');

clustattrib.nodata = 1;
if (resetaxes)
    plotgraph(handles);
end
%---------------------------------------------
function saveasMenu_Callback(hObject,handles)

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;

%o = analyzeoverlap(0);
% if (~isempty(o))
%    warndlg('You have overlapping clusters!','Warning');
% end
if (~isempty(clustattrib.currentfilename))
    suggestname = clustattrib.currentfilename;
else
    suggestname = ['matclust_',clustattrib.currentparamfilename];
end
if (ispc)
    [filename, pathname] = uiputfile('*.mat','Save MatClust file as...',suggestname);
else
    [filename, pathname] = filebrowse('save','filter','*.mat','title','Save MatClust file as...','suggestion',suggestname);
    
end
if (filename)
	currdir = pwd;
	
	cd(pathname);
	
	clustattrib.currentfilepath = pathname;
	clustattrib.currentfilename = filename; %stores the name of the current MatClust file
	set(handles.figure1,'Name',['MATCLUST  ',pathname,filename]);
	set(figattrib.handles.statusText,'String','Saving...');
    set(figattrib.handles.statusText,'UserData','Saving...');
    drawnow;
	
    currentStates = clustattrib.states;
    
    %save with undo states removed to reduce filesize
    clustattrib.states = {};                      
	
    save(clustattrib.currentfilename,'clustdata','graphattrib','clustattrib','-v6');   
     
    %put undo states back in
    clustattrib.states = currentStates;
    clustattrib.newchanges = 0;
    
    tmploc = strfind(filename,'.mat')-1;
    if ~isempty(tmploc)
        filename = filename(1:tmploc);
    end
    %cd([figattrib.datafoldername,'M_dataopen']);
    cd(figattrib.datafoldername);
    cd('M_dataopen');
    if (ispc)
        system(['rename ',figattrib.openfiles{figattrib.currentopenfile},' ',filename]);
    elseif (isunix)
        system(['mv "',figattrib.openfiles{figattrib.currentopenfile},'" "',filename,'"']);
    end
    figattrib.openfiles{figattrib.currentopenfile} = filename;
        
    cd(currdir);
    updatefileswitcher(handles);
    clustattrib.newchanges = 0;
    if ~isempty(clustattrib.currentfilename)
        set(handles.figure1,'Name',['MATCLUST  ',clustattrib.currentfilepath,clustattrib.currentfilename]);
    else
        set(handles.figure1,'Name',['MATCLUST  ',figattrib.openfiles{figattrib.currentopenfile}]);
    end
    
    set(figattrib.handles.statusText,'String','Ready');
    set(figattrib.handles.statusText,'UserData','Ready');
    addnewstate('SAVER',handles);
    
end
%----------------------------------------------
function saveMenu_Callback(hObject, handles)

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;

if (~isempty(clustattrib.currentfilename))
%    o = analyzeoverlap(0);
%     if (~isempty(o))
%         warndlg('You have overlapping clusters!','Warning');
%     end

    set(figattrib.handles.statusText,'String','Saving...');
    set(figattrib.handles.statusText,'UserData','Saving...');
    drawnow;
    cd(clustattrib.currentfilepath);
    
    
    currentStates = clustattrib.states;
    
    %save with undo states removed to reduce filesize
    clustattrib.states = {};                      
	save(clustattrib.currentfilename,'clustdata','graphattrib','clustattrib'); 
    
    %put undo states back in
    clustattrib.states = currentStates;
    clustattrib.newchanges = 0;
    
    if ~isempty(clustattrib.currentfilename)
        set(handles.figure1,'Name',['MATCLUST  ',clustattrib.currentfilepath,clustattrib.currentfilename]);
    else
        set(handles.figure1,'Name',['MATCLUST  ',figattrib.openfiles{figattrib.currentopenfile}]);
    end
    set(figattrib.handles.statusText,'String','Ready');
    set(figattrib.handles.statusText,'UserData','Ready');
   
    addnewstate('SAVER',handles);
else
    saveasMenu_Callback(hObject,handles)
end

%------------------------------------------------
function openclustfile_Callback(hObject,handles)

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;


if (ispc)
    [filename, pathname] = uigetfile('*.mat','Open MatClust file');
else
    [filename, pathname] = filebrowse('open','filter','*.mat','title','Open MatClust file');
end
if ((filename))	
    
    currdir = pwd;
	cd(pathname);
    fileinside = [];
    try
        fileinside = whos('-file', filename);
    end
    if ~((~isempty(fileinside)) & (length(fileinside)==3) & ((strcmp(fileinside(1).name,'clustdata'))|(strcmp(fileinside(1).name,'clustattrib'))|(strcmp(fileinside(1).name,'graphattrib'))))
        cd(currdir);
        error(['Not a valid MATCLUST file: ',filename]);
    end
    set(figattrib.handles.statusText,'String','Loading data...');
    set(figattrib.handles.statusText,'UserData','Loading data...');
    drawnow;
    if ~(clustattrib.nodata)
        savefigstate
    end
    clearworkspace(1,1);
    clearhistory;
    graphwindow = graphattrib.graphwindow;
	currclust = clustattrib.currclust;
	Himage = graphattrib.Himage;
	backgroundcolor = graphattrib.backgroundcolor;
    resolutionfactor = graphattrib.resolutionfactor;
    c0color = clustattrib.cluster0attrib.color;
    cluster0button = clustattrib.cluster0button;
    
	load(filename);
	
	graphattrib.Himage = Himage;
	clustattrib.currclust = currclust;
	graphattrib.graphwindow = graphwindow;
    graphattrib.backgroundcolor = backgroundcolor;
    graphattrib.resolutionfactor = resolutionfactor;
	clustattrib.cluster0attrib.color = c0color;
    clustattrib.cluster0button = cluster0button;
    clustattrib.currentfilepath = pathname;
    clustdata.timefiltersOn = zeros(32,1); %which time filters are currently on
    clustdata.timefiltersOn(1) = 1;
    clustdata.otherfiltersOn = zeros(32,1); %which 'other' filters are currently on
    clustdata.otherfiltersOn(1) = 1;
	if ~isfield(clustattrib,'dependencies')
        clustattrib.dependencies = false(size(figattrib.mixcolor,1));
    else
        if ((size(clustattrib.dependencies,1) < size(figattrib.mixcolor,1)) | (size(clustattrib.dependencies,2) < size(figattrib.mixcolor,1)))
            clustattrib.dependencies(size(figattrib.mixcolor,1),size(figattrib.mixcolor,1)) = 0;
        end
    end
	clustattrib.nodata = 0;
	open = OpenNewDataFolder;
    if (open)
        currdir = pwd;
        %cd([figattrib.datafoldername,'M_dataopen']);
        cd(figattrib.datafoldername);
        cd('M_dataopen');
        cd(figattrib.openfiles{figattrib.currentopenfile});
        load currstate;
        load states;
        clustattrib.states = states;
        maxstate = length(states);
        %load maxstate;
        cd(currdir);
        openstate(figattrib.openfiles{figattrib.currentopenfile},currstate,0,handles);
        if (currstate == 1)
            set(handles.undoMenu,'Label','Undo');
            set(handles.undoMenu,'Enable','off');
        else
            set(handles.undoMenu,'Label',['Undo ',clustattrib.lastaction]);
            set(handles.undoMenu,'Enable','on');
        end    
        if (currstate == maxstate)
            set(handles.redoMenu,'Enable','off');       
        else
            set(handles.redoMenu,'Enable','on');  
        end    
    else
        currdir = pwd;
        %cd([figattrib.datafoldername,'M_dataopen']);
        cd(figattrib.datafoldername);
        cd('M_dataopen');
        cd(figattrib.openfiles{figattrib.currentopenfile});
        currstate = 1;
        %maxstate = 1;
        paramdata = clustdata.params;
        save currstate currstate;
        %save maxstate maxstate;
        save paramdata paramdata;
        %clustdata = rmfield(clustdata,'params');
        baddim = find(diff(clustdata.datarange)==0);
        clustdata.datarange(1:2,baddim) = [clustdata.datarange(1,baddim)-1;clustdata.datarange(2,baddim)+1]; 
        clustattrib.newchanges = 1; %although no new changes have been made, this will alert new changes on the event of an undo
        
        clustattrib.states{1}.clustdata = rmfield(clustdata,'params');
        clustattrib.states{1}.clustattrib = rmfield(clustattrib,'states');
        clustattrib.states{1}.graphattrib = graphattrib;
        states = clustattrib.states;
        save states states;
        clustattrib.newchanges = 0;
        %clustdata.params = paramdata; 
        cd(currdir);            
        set(handles.undoMenu,'Label','Undo');
        set(handles.undoMenu,'Enable','off'); 
        set(handles.redoMenu,'Enable','off');
        set(handles.file_closeMenu,'Enable','on');
        set(handles.file_saveasMenu,'Enable','on');
        set(handles.file_saveMenu,'Enable','on');
        set(handles.file_curropenMenu,'Enable','on');
        set(handles.file_exportMenu,'Enable','on');
        set(handles.listbox1,'Enable','on');
        set(handles.listbox2,'Enable','on');
        set(handles.timefilterList,'Enable','on');
        set(handles.otherfilterList,'Enable','on');
        set(handles.paramaddContext1,'Enable','on');
        set(handles.paramaddContext2,'Enable','on');
        set(handles.parameditContext1,'Enable','on');
        set(handles.parameditContext2,'Enable','on');
        set(handles.paramdeleteContext1,'Enable','off');
        set(handles.paramdeleteContext2,'Enable','off');
        clustattrib.currstate = currstate;
    end
	
    fillworkspace(handles,1);
    updatefileswitcher(handles);
    plotgraph(handles);  
    setclustsizeinfo;
    if (ismember(clustattrib.currclust,clustattrib.clustersOn))
        set(handles.mainClustMenu,'Enable','on');
        for menuNum = 1:length(handles.mainClusterMenu)
            set(handles.mainClusterMenu(menuNum),'UserData',clustattrib.currclust)
        end
    else
        set(handles.mainClustMenu,'Enable','off');
    end
end

%-----------------------------------------------
function openstate(datafolder,statenum,filterpersist,handles)

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;



graphwindow = graphattrib.graphwindow;
currclust = clustattrib.currclust;
Himage = graphattrib.Himage;
backgroundcolor = graphattrib.backgroundcolor;
resolutionfactor = graphattrib.resolutionfactor;
c0color = clustattrib.cluster0attrib.color;
cluster0button = clustattrib.cluster0button;
timefilternames = clustdata.timefilternames;
otherfilternames = clustdata.otherfilternames;
timefiltersOn = clustdata.timefiltersOn;
otherfiltersOn = clustdata.otherfiltersOn;
hiddenclusters = clustattrib.hiddenclusters;
show = clustattrib.cluster0attrib.show;
currdir = pwd;

%cd([figattrib.datafoldername,'M_dataopen']);
cd(figattrib.datafoldername);
cd('M_dataopen');
cd(datafolder);

%load(['state',num2str(statenum)]);

load states;
load paramdata;

clustdata = states{statenum}.clustdata;
clustattrib = states{statenum}.clustattrib;
graphattrib = states{statenum}.graphattrib;
clustattrib.states = states;
clustdata.params = paramdata;

cd(currdir);
graphattrib.Himage = Himage;
clustattrib.currclust = currclust;
graphattrib.graphwindow = graphwindow;
graphattrib.backgroundcolor = backgroundcolor;
graphattrib.resolutionfactor = resolutionfactor;
clustattrib.cluster0attrib.color = c0color;
clustattrib.cluster0button = cluster0button;
if (filterpersist)
	clustdata.timefilternames = timefilternames;
	clustdata.otherfilternames = otherfilternames;
	clustdata.timefiltersOn = timefiltersOn;
	clustdata.otherfiltersOn = otherfiltersOn;
    clustattrib.hiddenclusters = hiddenclusters;
    clustattrib.cluster0attrib.show = show;
    FilterPoints
end
clustattrib.currstate = statenum;
clustattrib.nodata = 0;

%-----------------------------------------------
function fillworkspace(handles,resetaxeslists,KeepCurrentHiddenStates)

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;


if (nargin < 3)
    KeepCurrentHiddenStates = 0;
end
set(handles.timefilterList,'String',clustdata.timefilternames);
set(handles.otherfilterList,'String',clustdata.otherfilternames);
set(handles.timefilterList,'UserData',[1 0]);
set(handles.otherfilterList,'UserData',[1 0]);
set(handles.listbox1,'String',clustdata.names);
set(handles.listbox2,'String',clustdata.names);
if (resetaxeslists == 1) %reset the figure settings to startup values
    set(handles.listbox1,'Value',1);
    set(handles.listbox2,'Value',2);
    set(handles.listbox1,'UserData',1);
    set(handles.listbox2,'UserData',2);
    resetfilters(0);
    FilterPoints;
	hidevalue = 0;
	for i = 1:size(figattrib.clustcontrol,1)
         clusthideval(i,1) = {hidevalue};
	end
	set(figattrib.clustcontrol(:,3),{'Value'},clusthideval);
	clustattrib.hiddenclusters(:,1) = 0;
	excludevalue = 0;
	for i = 1:size(figattrib.clustcontrol,1)
         clustexcludeval(i,1) = {excludevalue};
	end
	set(figattrib.clustcontrol(:,4),{'Value'},clustexcludeval);
	clustattrib.hiddenclusters(:,2) = 0;  
    clustattrib.cluster0attrib.show = 1;
    
elseif (resetaxeslists == 0) %use figure settings for the file that is being opened
     currdir = pwd;
 	%cd([figattrib.datafoldername,'M_dataopen']);
 	cd(figattrib.datafoldername);
    cd('M_dataopen');
    cd(figattrib.openfiles{figattrib.currentopenfile});
 	load figstate;
 	
 	set(handles.listbox1,'Value', x);
 	set(handles.listbox1,'UserData', x);
 	set(handles.listbox2,'Value',y);
 	set(handles.listbox2,'UserData',y);
 	clustattrib.hiddenclusters = hiddenclusters;
 	for i = 1:size(clustattrib.hiddenclusters,1)
         set(figattrib.clustcontrol(i,3),'Value',hiddenclusters(i,1)); 
         set(figattrib.clustcontrol(i,4),'Value',hiddenclusters(i,2)); 
 	end
 	clustattrib.cluster0attrib.show = show;
 	clustdata.timefiltersOn = timefiltersOn;
 	clustdata.otherfiltersOn = otherfiltersOn;
    resetfilters(1);
    FilterPoints;
    cd(currdir);
elseif (resetaxeslists == 2) %dont change figure settings
    
end

if (clustattrib.cluster0attrib.show)
    set(handles.cluster0button,'String','Cluster 0: visible');	
    set(handles.cluster0button,'Value',1);	
else
    set(handles.cluster0button,'String','Cluster 0: hidden');	
    set(handles.cluster0button,'Value',0);
end
    
if ~isempty(clustattrib.currentfilename)
    if (clustattrib.newchanges)
        set(handles.figure1,'Name',['MATCLUST  ',clustattrib.currentfilepath,clustattrib.currentfilename,'*']);
    else
        set(handles.figure1,'Name',['MATCLUST  ',clustattrib.currentfilepath,clustattrib.currentfilename]);
    end    
else
    
    if (clustattrib.newchanges)
        set(handles.figure1,'Name',['MATCLUST  ',figattrib.openfiles{figattrib.currentopenfile},'*']);
    else
        set(handles.figure1,'Name',['MATCLUST  ',figattrib.openfiles{figattrib.currentopenfile}]);
    end        
end

axes(graphattrib.graphwindow);   

for clustnum = clustattrib.clustersOn' 

   
    turnoncluster(clustnum,handles,KeepCurrentHiddenStates)
    
    def = clustattrib.clusters{clustnum}.defineaxes; 
    for i = 1:size(def,1)
        
        a1 = def(i,1);
        a2 = def(i,2);
        findex = def(i,3);
        
        graphattrib.polyg(a1,a2,findex).lines{clustnum} = line(graphattrib.polyg(a1,a2,findex).vertices{clustnum}(:,1),graphattrib.polyg(a1,a2,findex).vertices{clustnum}(:,2));
        
        set(graphattrib.polyg(a1,a2,findex).lines{clustnum},'ButtonDownFcn','matclust(''polygon_ButtonDownFcn'',gcbo,guidata(gcbo))');
        
        set(graphattrib.polyg(a1,a2,findex).lines{clustnum},'UserData',[a1 a2 findex clustnum]);
        set(graphattrib.polyg(a1,a2,findex).lines{clustnum},'Color',figattrib.mixcolor(clustnum,:));
        set(graphattrib.polyg(a1,a2,findex).lines{clustnum},'UIContextMenu',handles.polycontext);
        
        if (graphattrib.polyg(a1,a2,findex).type{clustnum} == 2)                   
            set(graphattrib.polyg(a1,a2,findex).lines{clustnum},'LineStyle','-.','LineWidth',2); 
        end
        
        curra1 = get(handles.listbox1,'Value');
		curra2 = get(handles.listbox2,'Value');
		[currfindex, filterline] = findfilterindex(handles);
		%hide the currently displayed polygons
		
        if ~isequal([curra1 curra2 currfindex],[a1 a2 findex])
            set(graphattrib.polyg(a1,a2,findex).lines{clustnum},'Visible','off');
        end
        
    end
    
end

updateclustinfo(handles);
if (ismember(clustattrib.currclust,clustattrib.clustersOn))
    set(handles.mainClustMenu,'Enable','on');
    for menuNum = 1:length(handles.mainClusterMenu)
        set(handles.mainClusterMenu(menuNum),'UserData',clustattrib.currclust)
    end
else
    set(handles.mainClustMenu,'Enable','off');
end
%-------------------------------------------------
function alreadyopen = OpenNewDataFolder

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;

currdir = pwd;
%cd([figattrib.datafoldername,'M_dataopen']);
cd(figattrib.datafoldername);
cd('M_dataopen');
dirnames = dir;
figattrib.openfiles = [];

currnamecount = 1;
untitled = 0;
alreadyopen = 0;
filenum = [];
count = 1;

if ~isempty(clustattrib.currentfilename)
    tmploc = strfind(clustattrib.currentfilename,'.mat')-1;
    if ~isempty(tmploc)
        datafolder = clustattrib.currentfilename(1:tmploc);
    else
        datafolder = clustattrib.currentfilename;
    end
else
    datafolder = 'Untitled';    
    untitled = 1;
end

if (length(dirnames)>2)
    for i = 3:length(dirnames)
        if (dirnames(i).isdir)
            tmpname = dirnames(i).name;
            figattrib.openfiles{count} = tmpname;
            if (untitled)
                tmploc = strfind(tmpname,datafolder);
                if ~isempty(tmploc)
                    currnamecount = currnamecount+1;
                end
            else
                if ~(alreadyopen)
                    alreadyopen = strcmp(datafolder,tmpname);
                    filenum = count;
                end
            end        
            count = count+1;    
        end
    end    
end

if (untitled)
    datafolder = ['Untitled-',num2str(currnamecount)];
    figattrib.openfiles{end+1} = datafolder;
    figattrib.currentopenfile = length(figattrib.openfiles);
    mkdir(datafolder);
else
    if (alreadyopen)
        figattrib.currentopenfile = filenum;
    else
       figattrib.openfiles{end+1} = datafolder;
       figattrib.currentopenfile = length(figattrib.openfiles);
       mkdir(datafolder); 
   end
end     

cd(currdir);
%----------------------------------------------------
function changefilefocus(foldername,foldernum,resetaxes,resetfilt,filterpersist,handles)

%used to switch between currently open files. 
%foldername is the folder in M_dataopen where the to-be opened file is
%stored. resetaxes, resetfilt, filterpersist are either 1 or 0 describing
%how the axes and filter lists will behave when the switch is made

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;


savefigstate
currdir = pwd;
%cd([figattrib.datafoldername,'M_dataopen']);
cd(figattrib.datafoldername);
cd('M_dataopen');
cd(foldername);
load currstate;

cd(currdir);
figattrib.currentopenfile = foldernum;
for i = 1:length(figattrib.switchfileMenu)
     set(figattrib.switchfileMenu(i),'Checked','off');
end

set(figattrib.switchfileMenu(foldernum),'Checked','on');

clearworkspace(resetaxes,resetfilt);
clearhistory;
openstate(foldername,currstate,filterpersist,handles);
fillworkspace(handles,resetaxes);
if (currstate == 1)
    set(handles.undoMenu,'Label','Undo');
    set(handles.undoMenu,'Enable','off');
else
    set(handles.undoMenu,'Label',['Undo ',clustattrib.lastaction]);
    set(handles.undoMenu,'Enable','on');
end    
if (currstate == length(clustattrib.states))
    set(handles.redoMenu,'Enable','off');       
else
    set(handles.redoMenu,'Enable','on');  
end   

cd(currdir);
plotgraph(handles);   
%------------------------------------------------
function changefilefocus_after_close(foldername,foldernum,resetaxes,resetfilt,filterpersist,handles)

%used to switc to another open file after a file is closed. 
%foldername is the folder in M_dataopen where the to-be opened file is
%stored. resetaxes, resetfilt, filterpersist are either 1 or 0 describing
%how the axes and filter lists will behave when the switch is made

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;



currdir = pwd;
%cd([figattrib.datafoldername,'M_dataopen']);
cd(figattrib.datafoldername);
cd('M_dataopen');
cd(foldername);
load currstate;

cd(currdir);
figattrib.currentopenfile = foldernum;
for i = 1:length(figattrib.switchfileMenu)
     set(figattrib.switchfileMenu(i),'Checked','off');
end

set(figattrib.switchfileMenu(foldernum),'Checked','on');

clearworkspace(resetaxes,resetfilt);
clearhistory;
openstate(foldername,currstate,filterpersist,handles);
fillworkspace(handles,resetaxes);
if (currstate == 1)
    set(handles.undoMenu,'Label','Undo');
    set(handles.undoMenu,'Enable','off');
else
    set(handles.undoMenu,'Label',['Undo ',clustattrib.lastaction]);
    set(handles.undoMenu,'Enable','on');
end    
if (currstate == length(clustattrib.states))
    set(handles.redoMenu,'Enable','off');       
else
    set(handles.redoMenu,'Enable','on');  
end   

cd(currdir);
plotgraph(handles);   
%---------------------------------------
function updatefileswitcher(handles)

%updates the 'currently open files' menu

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;


if ~isempty(figattrib.switchfileMenu)
    for i = 1:length(figattrib.switchfileMenu)
        delete(figattrib.switchfileMenu(i));
    end
end
figattrib.switchfileMenu = [];


if ~isempty(figattrib.openfiles)
    for i = 1:length(figattrib.openfiles)
        figattrib.switchfileMenu(i) = uimenu(handles.file_curropenMenu,'Label',figattrib.openfiles{i}, ... 
            'Callback',['matclust(''changefilefocus'',''',figattrib.openfiles{i},''',',num2str(i),',0,1,0,guidata(gcbo))']);
        if (figattrib.currentopenfile == i)
            set(figattrib.switchfileMenu(i),'Checked','on');
        end
    end
end
%---------------------------------------
function addnewstate(action,handles)

%this function adds a new .mat file in the currently open folder in
%M-dataopen.  These files are titled state[statenum], where the highest
%number is the most recent.  When the number of states has reached the
%allowed maximum, the files are renamed so that the so that the oldest file
%is still numbered 1, and any older files are deleted.  Each state also has
%the action perforemed saved in clustattrib.lastaction

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;


if (~(strcmp(action,'zoom')) & ~(strcmp(action,'travel')))
    
    MaxAllowedState = figattrib.maxundos;
    currdir = pwd;
    %cd([figattrib.datafoldername,'M_dataopen']);
    cd(figattrib.datafoldername);
    cd('M_dataopen');
    cd(figattrib.openfiles{figattrib.currentopenfile});

     if isempty(strfind(action, 'SAVER'))
        clustattrib.lastaction = action;
        clustattrib.newchanges = 1;
        if ~isempty(clustattrib.currentfilename)
            set(handles.figure1,'Name',['MATCLUST  ',clustattrib.currentfilepath,clustattrib.currentfilename,'*']);
        else
            set(handles.figure1,'Name',['MATCLUST  ',figattrib.openfiles{figattrib.currentopenfile},'*']);
        end
        takeoff = 0;
        if (clustattrib.currstate >= MaxAllowedState)
           takeoff = clustattrib.currstate - MaxAllowedState + 1;
        end
        clustattrib.states = clustattrib.states(takeoff+1:clustattrib.currstate);
        clustattrib.currstate = length(clustattrib.states) + 1;
        set(handles.undoMenu,'Label',['Undo ',clustattrib.lastaction]);
        set(handles.redoMenu,'Enable','off');
        set(handles.undoMenu,'Enable','on');  
        if ~isempty(strfind(clustattrib.lastaction, 'parameter'))
            paramdata = clustdata.params;
            save paramdata paramdata;
            set(handles.undoMenu,'Label','Undo');
            set(handles.undoMenu,'Enable','off');  
        end
     
     else
         %user just saved- just append info to the current state
         
     end
    

   
        
    
    clustattrib.states{clustattrib.currstate}.clustattrib = rmfield(clustattrib,'states');
    clustattrib.states{clustattrib.currstate}.clustdata = rmfield(clustdata,'params');
    clustattrib.states{clustattrib.currstate}.graphattrib = graphattrib;

    
    cd(currdir);
    
else
    clustattrib.states{clustattrib.currstate}.graphattrib = graphattrib;

end
    
%-------------------------------------------------
function undoMenu_Callback(hObject,handles)

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;

resetaxes = 0;
resetaxes2 = 2;
resetfilt = 0;
filterpersist = 1;

currfilename = clustattrib.currentfilename;
currfilepath = clustattrib.currentfilepath;
set(figattrib.handles.statusText,'String','Loading data...');
set(figattrib.handles.statusText,'UserData','Loading data...');

timefiltersOn = clustdata.timefiltersOn;
otherfiltersOn = clustdata.otherfiltersOn;
if ~isempty(strfind(clustattrib.lastaction, 'filter'))
    resetfilt = 1;
    filterpersist = 0;
end


viewbox = graphattrib.viewbox;
params = clustdata.params;
states = clustattrib.states;
statenum = clustattrib.currstate;
hiddenclusters = clustattrib.hiddenclusters;
clearworkspace(resetaxes,resetfilt);

graphwindow = graphattrib.graphwindow;
currclust = clustattrib.currclust;
Himage = graphattrib.Himage;
backgroundcolor = graphattrib.backgroundcolor;
resolutionfactor = graphattrib.resolutionfactor;
c0color = clustattrib.cluster0attrib.color;
cluster0button = clustattrib.cluster0button;
timefilternames = clustdata.timefilternames;
otherfilternames = clustdata.otherfilternames;


show = clustattrib.cluster0attrib.show;

clustdata = states{statenum-1}.clustdata;
clustattrib = states{statenum-1}.clustattrib;
graphattrib = states{statenum-1}.graphattrib;
clustattrib.states = states;
clustdata.params = params;

clustattrib.currclust = currclust;
graphattrib.backgroundcolor = backgroundcolor;
graphattrib.resolutionfactor = resolutionfactor;
clustattrib.cluster0attrib.color = c0color;
clustattrib.cluster0button = cluster0button;
if (filterpersist)
	clustdata.timefilternames = timefilternames;
	clustdata.otherfilternames = otherfilternames;
	clustdata.timefiltersOn = timefiltersOn;
	clustdata.otherfiltersOn = otherfiltersOn;
    clustattrib.hiddenclusters = hiddenclusters;
    clustattrib.cluster0attrib.show = show;
    FilterPoints
end

clustattrib.currstate = statenum - 1;
clustattrib.nodata = 0;
%openstate(foldername,currstate,filterpersist,handles);

if ~isempty(currfilename)
    clustattrib.currentfilename = currfilename;
    clustattrib.states{clustattrib.currstate}.clustattrib.currentfilename = currfilename;
    clustattrib.currentfilepath = currfilepath;
end


fillworkspace(handles,resetaxes2,1);

if (~filterpersist)
    clustdata.timefiltersOn = timefiltersOn;
    clustdata.otherfiltersOn = otherfiltersOn;
    clustattrib.hiddenclusters = hiddenclusters;
    clustattrib.cluster0attrib.show = show;
    
    resetfilters(2)
    FilterPoints
end
  
graphattrib.viewbox = viewbox;
clustattrib.states{clustattrib.currstate}.graphattrib.viewbox = viewbox;

plotgraph(handles);   
set(handles.redoMenu,'Enable','on');
if (clustattrib.currstate == 1)
    set(handles.undoMenu,'Label','Undo');
    set(handles.undoMenu,'Enable','off');
else
    set(handles.undoMenu,'Label',['Undo ',clustattrib.lastaction]);
end
if ~isempty(strfind(clustattrib.lastaction, 'parameter'))
    set(handles.undoMenu,'Label','Undo');
    set(handles.undoMenu,'Enable','off');  
end


%----------------------------------------------
function redoMenu_Callback(hObject,handles)

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;

set(figattrib.handles.statusText,'String','Loading data...');
set(figattrib.handles.statusText,'UserData','Loading data...');
currfilename = clustattrib.currentfilename;
currfilepath = clustattrib.currentfilepath;
resetaxes = 0;
resetaxes2 = 2;
resetfilt = 0;
filterpersist = 1;
viewbox = graphattrib.viewbox;
timefiltersOn = clustdata.timefiltersOn;
otherfiltersOn = clustdata.otherfiltersOn;



params = clustdata.params;
states = clustattrib.states;
statenum = clustattrib.currstate;
nextaction = states{statenum+1}.clustattrib.lastaction;
if ~isempty(strfind(nextaction, 'filter'))
    resetfilt = 1;
    filterpersist = 0;
end
hiddenclusters = clustattrib.hiddenclusters;
clearworkspace(resetaxes,resetfilt);

graphwindow = graphattrib.graphwindow;
currclust = clustattrib.currclust;
Himage = graphattrib.Himage;
backgroundcolor = graphattrib.backgroundcolor;
resolutionfactor = graphattrib.resolutionfactor;
c0color = clustattrib.cluster0attrib.color;
cluster0button = clustattrib.cluster0button;
timefilternames = clustdata.timefilternames;
otherfilternames = clustdata.otherfilternames;
show = clustattrib.cluster0attrib.show;

clustdata = states{statenum+1}.clustdata;
clustattrib = states{statenum+1}.clustattrib;
graphattrib = states{statenum+1}.graphattrib;
clustattrib.states = states;
clustdata.params = params;

clustattrib.currclust = currclust;
graphattrib.backgroundcolor = backgroundcolor;
graphattrib.resolutionfactor = resolutionfactor;
clustattrib.cluster0attrib.color = c0color;
clustattrib.cluster0button = cluster0button;
if (filterpersist)
	clustdata.timefilternames = timefilternames;
	clustdata.otherfilternames = otherfilternames;
	clustdata.timefiltersOn = timefiltersOn;
	clustdata.otherfiltersOn = otherfiltersOn;
    clustattrib.hiddenclusters = hiddenclusters;
    clustattrib.cluster0attrib.show = show;
    FilterPoints
end
clustattrib.currstate = statenum + 1;
clustattrib.nodata = 0;
%openstate(foldername,currstate,filterpersist,handles);
if ~isempty(currfilename)
    clustattrib.currentfilename = currfilename;
    clustattrib.states{clustattrib.currstate}.clustattrib.currentfilename = currfilename;
    clustattrib.currentfilepath = currfilepath;
end
fillworkspace(handles,resetaxes2,1);
if (~filterpersist)
    clustdata.timefiltersOn = timefiltersOn;
    clustdata.otherfiltersOn = otherfiltersOn;
    clustattrib.hiddenclusters = hiddenclusters;
    clustattrib.cluster0attrib.show = show;
    
    resetfilters(2)
    FilterPoints
end
graphattrib.viewbox = viewbox;
clustattrib.states{clustattrib.currstate}.graphattrib.viewbox = viewbox;

plotgraph(handles);   

if (clustattrib.currstate == length(states))
    set(handles.redoMenu,'Enable','off');
end

set(handles.undoMenu,'Label',['Undo ',clustattrib.lastaction]);
set(handles.undoMenu,'Enable','on');  

%----------------------------------------
function savefigstate()

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;

handles = figattrib.handles;
currdir = pwd;
%cd([figattrib.datafoldername,'M_dataopen']);
cd(figattrib.datafoldername);
cd('M_dataopen');
cd(figattrib.openfiles{figattrib.currentopenfile});

x = get(handles.listbox1,'Value');
y = get(handles.listbox2,'Value');
hiddenclusters = clustattrib.hiddenclusters;
show = clustattrib.cluster0attrib.show;
timefiltersOn = clustdata.timefiltersOn;
otherfiltersOn = clustdata.otherfiltersOn;


save figstate x y hiddenclusters show timefiltersOn otherfiltersOn;

states = clustattrib.states;
currstate = clustattrib.currstate;
save states states;
save currstate currstate; 
cd(currdir);

%----------------------------------------------
function exporttimefiltMenu_Callback(handles)

global clustdata

currdir = pwd;
ranges = clustdata.timefilterranges;
names = clustdata.timefilternames;
cd ..
if (ispc)
    [filename, pathname] = uiputfile('*.mat','Export time filters file as...','times');
else
    [filename, pathname] = filebrowse('save','filter','*.mat','title','Export time filters file as...','suggestion','times.mat');
end
if (filename)
	%currdir = pwd;
	cd(pathname);
    save(filename,'ranges','names');
    %cd(currdir)
end
cd(currdir)
%-----------------------------------------------
function importtimefiltMenu_Callback(handles)

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;

%if(isempty(clustattrib.clustersOn))
if (sum((clustdata.timefilterranges(:,1)~=0)|(clustdata.timefilterranges(:,2)~=0)) < 2)
   
	if (ispc)
        [filename, pathname] = uigetfile('*.mat','Find time filter file');
    else
        [filename, pathname] = filebrowse('open','filter','*.mat','title','Find time filter file');
    end
	
	if (filename)        
        importtimefilt(pathname, filename);               
	end
else
    msgbox('This operation can only be done if no time filters are currently defined.  Delete all existing time filters.');
    %beep
    %disp('This operation can only be done if no time filters are currently defined')
end
%------------------------------------------
function importtimefilt(pathname, filename)

% imports a time filters file in the pathname directory

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;

if(clustattrib.nodata)
     error('You must have data loaded before importing a times filter file');
end
handles = figattrib.handles;
currdir = pwd;
try
    cd(pathname);
    a = whos('-file',filename);
    if ~(strcmp(a(1).name,'ranges') | strcmp(a(2).name,'ranges'))
        error('---');
    end
    load(filename);
    cd(currdir);
    allpointsrange = clustdata.timefilterranges(1,:);
    ranges(1,1:2) = allpointsrange; 
    clustdata.timefilterranges = ranges;
    clustdata.timefilternames = names;
catch
    cd(currdir);
    error('Not a valid times filter file.');
end

set(handles.timefilterList,'String',names);
set(handles.timefilterList,'UserData',[1 0]);

clustdata.timefilters = int32(zeros(size(clustdata.params,1),1)); %contains the filtere points of up to 32 time filters
clustdata.timefilters = fastbitset(clustdata.timefilters,1,logical(1)); %the first filter is all times
clustdata.timefiltermemmap = zeros(32,1); %because the flter order can be changed by the user, this map translates the filter numbers to those stored in memory
clustdata.timefiltermemmap(1) = 1;
clustdata.timefiltersOn = clustdata.timefiltermemmap'; %which time filter is currently on
clustdata.timefiltermemmap(1) = 0;
for i = 1:size(ranges,1)

    passfilter = ((clustdata.params(:,1) >= ranges(i,1))&(clustdata.params(:,1) <= ranges(i,2)));
    if sum(passfilter)
        firstzero = min(find(clustdata.timefiltermemmap==0));
        clustdata.timefiltermemmap(firstzero) = i;
        clustdata.timefilters = fastbitset(clustdata.timefilters,firstzero,passfilter);
    end
end
%set(fighandles.mainfighandles.timefilterList,'Value',[]);
cleartimefilter;
clustdata.filtermemmap = [clustdata.timefiltermemmap; clustdata.otherfiltermemmap];

resetfilters(0);
addnewstate('import time filters',handles);
%------------------------------------------------------------------

function exportimageMenu_Callback(handles)
% function used to copy the current cluster image to a new figure

global graphattrib;
global clustattrib;
global figattrib;

%if (ispc)
    [filename, pathname] = uiputfile('*.png','Save image as...');
%else
%    [filename, pathname] = filebrowse('save','filter','*.bmp','title','Save image as...');   
%end


%newframe = getframe(gca);
 XT = get(gca,'XTick');
 YT = get(gca,'YTick');
 XTL = get(gca,'XTickLabel');
 YTL = flipud(get(gca,'YTickLabel'));
 currpos = get(gcf,'position');
%newimage = frame2im(newframe);
% tmpfig = figure;
% newaxes = axes;
% newimagehandle = image(newimage);


%newimage = flipud(get(graphattrib.Himage,'CData'));
newfig = figure;
newaxes = axes;
%set(newaxes,'YDir','reverse');
set(newaxes,'XTick',XT);
set(newaxes,'YTick',YT);
set(newaxes,'XTickLabel',XTL);
set(newaxes,'YTickLabel',flipud(YTL));
new_handle = copyobj(graphattrib.Himage,newaxes);
colormap([graphattrib.backgroundcolor;clustattrib.cluster0attrib.color;figattrib.mixcolor]);
set(newfig,'position',currpos);

 currdir = pwd;
 cd(pathname)
 print(['-f',num2str(gcf)],'-dpng',filename);
% imwrite(newimage,[graphattrib.backgroundcolor;clustattrib.cluster0attrib.color;figattrib.mixcolor],filename,'BMP')
 cd(currdir);
 delete(newfig);
%-----------------------------------------
function setclustsizeinfo()

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;

for i = 1:length(figattrib.mixcolor)
    if ismember(i,clustattrib.clustersOn)
        set(figattrib.clustcontrol(i,2),'TooltipString',['Cluster ',num2str(i),': ',num2str(length(clustattrib.clusters{i}.index)),' points']);
    else
        set(figattrib.clustcontrol(i,2),'TooltipString',['Cluster ',num2str(i)]);
    end
end
%------------------------------------------
function exportdataMenu_Callback(handles)

global clustdata;
global graphattrib;
global clustattrib;
global figattrib;

o = analyzeoverlap(0);
if (~isempty(o))
    msgbox('You have overlapping clusters! Export aborted.','Error');
    return;
end
name = figattrib.openfiles{figattrib.currentopenfile};
clustparam = zeros(size(clustdata.params,1),1);
clustparam(:,2) = 1;

for i = 1:length(clustattrib.clusters)
    try
        clustparam(clustattrib.clusters{i}.index,1) = i;
    end
end

for i = 1:size(clustdata.timefilterranges,1)
    
    if (i >1)
        clustparam(find((clustdata.params(:,1)>=clustdata.timefilterranges(i,1))&(clustdata.params(:,1)<=clustdata.timefilterranges(i,2))),2) = i;
    end
end

save(['Export_',name],'clustparam');
%--------------------------------------------
function ExecuteTool(funcname);

global figattrib;

currdir = pwd;

cd(figattrib.foldername);
cd('Tools');
eval(['F = @',funcname,';']);
cd(currdir);
%eval(funcname);
feval(F);

%------------------------------------------
function ExecuteClustTool(funcname,clustnum);

global figattrib;

currdir = pwd;

cd(figattrib.foldername);
cd('ClustTools');
eval(['F = @',funcname,';']);
cd(currdir);
%eval([funcname,'(',num2str(clustnum),');']);
feval(F,clustnum);

%------------------------------------------
function TFiltBehaveFcn(hObject, handles)
% controls the time filter behavior (a menu function)
global figattrib;

caller = get(hObject,'UserData');  %1 for 'one only' and 2 for 'allow multiple'
switch caller
    case 1
        set(handles.oneTFiltMenu,'Checked','On');
        set(handles.multTFiltMenu,'Checked','Off');
        set(handles.timeFiltOptButton,'Value',1);
        set(handles.timeFiltOptButton,'TooltipString','Time selection mode: one at a time');
        figattrib.tFiltAllowMultiple = 0;
    case 2
        set(handles.oneTFiltMenu,'Checked','Off');
        set(handles.multTFiltMenu,'Checked','On');
        set(handles.timeFiltOptButton,'Value',0);
        set(handles.timeFiltOptButton,'TooltipString','Time selection mode: allow multiple');
        figattrib.tFiltAllowMultiple = 1;
    case 3 %the mode button was used instead of the menu option
        mode = get(hObject,'Value');
        if ~(mode) %mode button is not on (multiple selection mode)
            set(handles.oneTFiltMenu,'Checked','Off');
            set(handles.multTFiltMenu,'Checked','On');
            set(handles.timeFiltOptButton,'TooltipString','Time selection mode: allow multiple');
            figattrib.tFiltAllowMultiple = 1;
        else %single selection mode
            set(handles.timeFiltOptButton,'TooltipString','Time selection mode: one at a time');
            set(handles.oneTFiltMenu,'Checked','On');
            set(handles.multTFiltMenu,'Checked','Off');
            figattrib.tFiltAllowMultiple = 0;
        end
                    
end
%-------------------------------------------
function TUnitsMenuFcn(hObject, handles)

global clustdata
global figattrib

answer = inputdlg('One second =','Time Units',1,{num2str(clustdata.UnitsPerSec)});
if (isempty(answer))
    return
end
UnitsPerSec = str2num(answer{1});
if (isempty(UnitsPerSec))
    return
end
clustdata.UnitsPerSec = UnitsPerSec;
currdir = pwd;
cd(figattrib.foldername);
load matclust_defaults;
matclust_defaults.UnitsPerSec = clustdata.UnitsPerSec;
save('matclust_defaults','matclust_defaults');
cd(currdir);
plotgraph(handles);
%------------------------------------------
function matclusthelpMenuFcn(handles)

global figattrib
currdir = pwd;
cd(figattrib.foldername);

open('matclusthelp.html');
cd(currdir);

%-------------------------------------------
function hidecurrentpolygons(handles)

global graphattrib;

a1 = get(handles.listbox1,'UserData');
a2 = get(handles.listbox2,'UserData');
[findex, filterline] = findfilterindex(handles);
%hide the currently displayed polygons
for fNum = 1:length(findex)
    try
        for i = 1:length(graphattrib.polyg(a1,a2,findex(fNum)).lines)
            try
                set(graphattrib.polyg(a1,a2,findex(fNum)).lines{i},'Visible','off');
                delete(graphattrib.polyg(a1,a2,findex(fNum)).highlight{i});    
            end
        end
    end
end

%---------------------------------------------
function showcurrentpolygons(handles)

global graphattrib;

a1 = get(handles.listbox1,'UserData');
a2 = get(handles.listbox2,'UserData');
[findex, filterline] = findfilterindex(handles);
%show the new polygons in in the active filter
for fNum = 1:length(findex)
    try
        for i = 1:length(graphattrib.polyg(a1,a2,findex(fNum)).lines)
            try
                set(graphattrib.polyg(a1,a2,findex(fNum)).lines{i},'Visible','on');        
            end
        end
    end    
end
%----------------------------------------------
function color = colorclip(color)

color(find(color < 0)) = 0;
%---------------------------------------------
function ReleaseFocus(fig)
%annoying work-around to give main fig focus

set(findobj(fig, 'Type', 'uicontrol'), 'Enable', 'off');
drawnow;
set(findobj(fig, 'Type', 'uicontrol'), 'Enable', 'on');
%----------------------------------------------
function highlightCurrentCluster(handles)

global graphattrib;
global clustattrib;

a1 = get(handles.listbox1,'UserData');
a2 = get(handles.listbox2,'UserData');
objecthandle = [];
[findex, filterline] = findfilterindex(handles);
for i = 1:length(findex)
    try
        objecthandle = [objecthandle graphattrib.polyg(a1,a2,findex(i)).lines{clustattrib.currclust}];
    end
end
set(objecthandle,'Marker','s'); %highlight the clicked polygon
graphattrib.currentpolyhighlight = objecthandle;

if (length(objecthandle) == 1)
    %set(handles.mainPolygonMenu,'Enable','on');
    set(handles.polydeletemenu,'Enable','on');
    set(handles.addPolygonFilterMainmenu,'Enable','on');
    set(handles.copyPolygonFilterMainmenu,'Enable','on');
    
    Child = get(handles.polycontext,'Children');
    for ch = 1:length(Child)
        set(Child(ch),'Enable','on');
    end
else
    %set(handles.mainPolygonMenu,'Enable','off');
    set(handles.polydeletemenu,'Enable','on');
    set(handles.addPolygonFilterMainmenu,'Enable','off');
    set(handles.copyPolygonFilterMainmenu,'Enable','off');
    Child = get(handles.polycontext,'Children');
    for ch = 1:length(Child)
        set(Child(ch),'Enable','off');
    end
end
%------------------------------------------------------
function resetButtonBackgrounds(handles)

global figattrib;

bcolor = figattrib.figureColor*255;
temppic = imread('poly','bmp');
temppic = changepicbackground(temppic,bcolor);
set(handles.polygonbutton,'CData',temppic);
temppic = imread('square','bmp');
temppic = changepicbackground(temppic,bcolor);
set(handles.squarebutton,'CData',temppic);
temppic = imread('arrow','bmp');
temppic = changepicbackground(temppic,bcolor);
set(handles.arrowbutton,'CData',temppic);
temppic = imread('mag','bmp');
temppic = changepicbackground(temppic,bcolor);
set(handles.magbutton,'CData',temppic);
temppic = imread('hand','bmp');
temppic = changepicbackground(temppic,bcolor);
set(handles.handbutton,'CData',temppic);
temppic = imread('wandpic','bmp');
temppic = changepicbackground(temppic,bcolor);
set(handles.wandbutton,'CData',temppic);

