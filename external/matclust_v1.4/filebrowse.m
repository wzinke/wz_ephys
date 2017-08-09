function [filename, path] = filebrowse(mode,varargin)
%[filename, path] = filebrowse(mode,properties)
%   
%   This function opens a file browsing dialog and allows the user to pick
%   a file.  There are two different modes set by the 'mode' input
%   variable.  The two modes are 'save' or 'open', and the dialog will
%   behave slightly different in the two modes.
%
%   There are also some properties which can be set:
%   'filter' - this sets the file filter (such as '*.m'). Default is '*'.
%   'title'  - sets the title of the dialog
%   'suggestion' - gives a suggested filename in the filename box of the dialog
%   'directory' - gives a starting directory
%

if (strcmp(lower(mode),'save')) %are we in save or open mode?
    filemode = 1;
    mode = 'Save';
elseif (strcmp(lower(mode),'open'))
    filemode = 2;
    mode = 'Open';
else
    error('Mode must be either ''save'' or ''open''');
end

filter = '*';
title = ['Choose name of file to ',lower(mode)];
suggestion = '';
basedir = pwd;
%this is where the optional inputs are set
for i = 1:2:length(varargin)
    switch (varargin{i})
        case 'filter'
            filter = varargin{i+1};
        case 'title'
            title = varargin{i+1};
        case 'suggestion'
            suggestion = varargin{i+1};
        case 'directory'
            basedir = varargin{i+1};
        otherwise
            error([varargin{i},' is not a valid option']);
    end
end
if (nargin < 2)
    filter = '*';
end
    
%create the figure
fighandle = figure('MenuBar','none','Color',get(0,'DefaultUicontrolBackgroundColor'),'Position',[400 400 500 300], ... 
        'Tag','eventviewer','NumberTitle','off','Resize','off','Name',title, ... 
        'CloseRequestFcn','filebrowse_controls(''CloseRequestFnc'',guidata(gcbo))', ...
        'WindowStyle','normal');
    
currdir = pwd;
try
    cd(basedir);
end

%find the files and directories to display in the listbox
diroutput = dir;
diroutput = diroutput(2:end);
dirlist = [];
count = 1;
for a = 1:length(diroutput)
    if diroutput(a).isdir
        dirlist{count} = ['<',diroutput(a).name,'>'];
        dirstruct(count) = diroutput(a);
        count = count+1;
    end
end
if (count > 1) %some folders are displayed
    dirstruct(count).isdir = -1;
    dirlist{count} = '==========================';
    count = count+1;
end
diroutput = dir(filter);
for a = 1:length(diroutput)  
    if ~diroutput(a).isdir     
        dirlist{count} = diroutput(a).name;
        dirstruct(count) = diroutput(a);
        count = count+1;
    end
end
    
%create the edit boxes   
e(1) = uicontrol('Tag','pathedit','Style','edit','Units','normalized','HorizontalAlignment','left','BackgroundColor',[1 1 1],'Position',[.02 .85 .96 .07],'String',pwd,'Callback','filebrowse_controls(''pathedit_Callback'',guidata(gcbo))');
e(2) = uicontrol('Tag','fileedit','Style','edit','Units','normalized','HorizontalAlignment','left','BackgroundColor',[1 1 1],'Position',[.15 .2 .83 .07],'String',suggestion,'Callback','filebrowse_controls(''fileedit_Callback'',guidata(gcbo))');
e(2) = uicontrol('Tag','filteredit','Style','edit','Units','normalized','HorizontalAlignment','left','BackgroundColor',[1 1 1],'Position',[.15 .1 .35 .07],'String',filter,'Callback','filebrowse_controls(''filteredit_Callback'',guidata(gcbo))');

%create text boxes
t(1) = uicontrol('Tag','text1','Style','text','Units','normalized','FontUnits','pixels','FontSize',11,'Position',[.02 .2 .1 .07],'String','Filename:','HorizontalAlignment','right');
t(2) = uicontrol('Tag','text2','Style','text','Units','normalized','FontUnits','pixels','FontSize',11,'Position',[.02 .1 .1 .07],'String','Filter:','HorizontalAlignment','right');

%create buttons
b(1) = uicontrol('Tag','okbutton','Style','pushbutton','Units','normalized','Position',[.88 .03 .10 .07],'String',mode,'Callback','filebrowse_controls(''okbutton_Callback'',guidata(gcbo))');
b(2) = uicontrol('Tag','cancelbutton','Style','pushbutton','Units','normalized','Position',[.775 .03 .10 .07],'String','Cancel','Callback','filebrowse_controls(''cancelbutton_Callback'',guidata(gcbo))');

%create the listbox
listcontrol(1) = uicontrol('Style','listbox','Tag','listbox1','BackgroundColor',[1 1 1],'ForegroundColor',[0 0 0], ...
'Units','normalized','Position',[.02 .3 .96 .5],'String',dirlist,'Max',1,'Min',1, ...
'Callback','filebrowse_controls(''listbox1_Callback'',guidata(gcbo))', 'Value',[1]);


handles = guihandles(fighandle);
handles.endvar = 0;
handles.filemode = filemode;
handles.filter = filter;
handles.dirlist = dirlist;
handles.dirstruct = dirstruct;
handles.path = pwd;
handles.filename = suggestion;
handles.fighandle = fighandle;
guidata(fighandle, handles);

%the gui will now force matlab to wait until one of the helper functions
%releases it
uiwait(fighandle);

%when the gui is released, we output the pathname and filename, and finally
%the figure is deleted
handles = guidata(fighandle);
path = handles.path;
if (isunix)
    if (~isempty(path) & ~(strcmp(path(end),'/')))
        path = [path,'/'];
    end
else
    if (~isempty(path) & ~(strcmp(path(end),'\')))
        path = [path,'\'];
    end
end

filename = handles.filename;
cd(currdir);
delete(fighandle);




