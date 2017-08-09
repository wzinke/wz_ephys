function setuplocal(varargin)

if (nargin == 0)
    start;
else
     feval(varargin{:});
     
end


function start

%mainfighandles = handles;

%mainfigpos = get(handles.figure1,'Position');

fighandle = figure('MenuBar','none','Color',get(0,'DefaultUicontrolBackgroundColor'),'Tag','LocalSetupFig','NumberTitle','off','Resize','off','Name','Matclust Local Setup', ... 
        'CloseRequestFcn','setuplocal(''LocalSetupFig_CloseRequestFnc'',gcbo,guidata(gcbo))', ...
        'Position', [400 400 240 200]);

e(1) = uicontrol('Tag','edit1','Style','edit','Position',[20 140 165 20],'BackgroundColor',[1 1 1],'Callback','setuplocal(''edit1_Callback'',gcbo,guidata(gcbo))');
e(2) = uicontrol('Tag','edit2','Style','edit','Position',[20 70 165 20],'BackgroundColor',[1 1 1],'Callback','setuplocal(''edit2_Callback'',gcbo,guidata(gcbo))');
%e(3) = uicontrol('Tag','edit3','Style','edit','Position',[70 52 115 20],'BackgroundColor',[1 1 1]);
b(1) = uicontrol('Tag','cancelbutton','Style','pushbutton','Position',[126 5 50 20],'String','Cancel','Callback','setuplocal(''cancelbutton_Callback'',gcbo,guidata(gcbo))');
b(2) = uicontrol('Tag','okbutton','Style','pushbutton','Position',[177 5 50 20],'String','Ok','Callback','setuplocal(''okbutton_Callback'',gcbo,guidata(gcbo))');
b(3) = uicontrol('Tag','lookbutton1','Style','pushbutton','Position',[187 140 25 20],'String','/','Callback','setuplocal(''lookbutton1_Callback'',gcbo,guidata(gcbo))');
b(4) = uicontrol('Tag','lookbutton2','Style','pushbutton','Position',[187 70 25 20],'String','/','Callback','setuplocal(''lookbutton2_Callback'',gcbo,guidata(gcbo))');

t(1) = uicontrol('Tag','text1','Style','text','Position',[15 160 200 30],'String','Where should the matclust folder be created?','HorizontalAlignment','left');
t(2) = uicontrol('Tag','text2','Style','text','Position',[15 90 200 30],'String','Give a directory on the local hardrive for temporary file writing','HorizontalAlignment','left');
%t(3) = uicontrol('Tag','text3','Style','text','Position',[15 50 50 20],'String','End time:','HorizontalAlignment','left');
%t(4) = uicontrol('Tag','text4','Style','text','Position',[5 180 200 20],'String','','HorizontalAlignment','left','ForegroundColor',[1 0 0],'BackgroundColor',[.925 .914 .847]);
fighandles = guihandles(fighandle);
fighandles.fighandle = fighandle;
%fighandles.mainfighandles = mainfighandles;
fighandles.pathname1 = [];
fighandles.pathname2 = [];
guidata(fighandle,fighandles);
%------------------------------------------------------------
function LocalSetupFig_CloseRequestFnc(hObject,fighandles)

delete(hObject);
%------------------------------------------------------------
function cancelbutton_Callback(hObject,fighandles)

delete(fighandles.LocalSetupFig);
%------------------------------------------------------------
function lookbutton1_Callback(hObject,fighandles)



dirname = uigetdir('','Where should the matclust folder be created?');

if ischar(dirname)
    set(fighandles.edit1,'String',dirname);
    fighandles.pathname1 = dirname;
end
guidata(fighandles.fighandle,fighandles);

%------------------------------------------------------------
function lookbutton2_Callback(hObject,fighandles)



dirname = uigetdir('','Give a directory on the local hardrive for temporary file writing');

if ischar(dirname)
    set(fighandles.edit2,'String',dirname);
    fighandles.pathname2 = dirname;
end
guidata(fighandles.fighandle,fighandles);
%---------------------------------------------------------

function okbutton_Callback(hObject,fighandles)

fighandles.pathname1
fighandles.pathname2

if (isempty(fighandles.pathname1)|isempty(fighandles.pathname2))
    'you must fill both fields'
    return    
end
 
currdir = pwd;

tmploc = strfind(which('setuplocal'),'.m')-11;
filelocation = which('setuplocal');
filelocation = filelocation(1:tmploc);  %contains the path to MatClust

cd(fighandles.pathname1);
mkdir matclust;

cd matclust;
matclustfolder = pwd;


mkdir ParameterPrograms;
cd ParameterPrograms;
paramdir = pwd;
cd ..
mkdir Tools;
cd Tools;
toolsdir = pwd;
cd ..
mkdir ClustTools;
cd ClustTools;
clustdir = pwd;
cd ..
mkdir Filters;
cd Filters;
filterdir = pwd;
cd ..

cd(filelocation);

copyfile('ParameterPrograms',paramdir);
copyfile('Tools',toolsdir);
copyfile('ClustTools',clustdir);
copyfile('Filters',filterdir);
copyfile('*.c',matclustfolder);

cd(filelocation)
load defaultsbin;
cd(fighandles.pathname1);
cd matclust;


matclust_defaults.Cluster0Color = defaults.Cluster0Color;
matclust_defaults.ClusterColors = defaults.ClusterColors;
matclust_defaults.GraphBackgroundColor = defaults.GraphBackgroundColor;
matclust_defaults.MaxUndos = defaults.MaxUndos;
matclust_defaults.ResFactor = defaults.ResFactor;
matclust_defaults.UnitsPerSec = defaults.UnitsPerSec;
matclust_defaults.DataFolder = fighandles.pathname2;
matclust_defaults.UserFolder = pwd;

save matclust_defaults matclust_defaults; 

cfiles = dir('*.c');
for i = 1:length(cfiles)
   mex(cfiles(i).name)
end
delete *.c;
cd(currdir);
delete(fighandles.LocalSetupFig);

%---------------------------------------
function edit1_Callback(hObject,fighandles)

fighandles.pathname1 = get(hObject,'String');
guidata(fighandles.fighandle,fighandles);
%---------------------------------------
function edit2_Callback(hObject,fighandles)

fighandles.pathname2 = get(hObject,'String');
guidata(fighandles.fighandle,fighandles);