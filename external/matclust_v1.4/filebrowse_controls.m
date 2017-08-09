function filebrowse_controls(varargin)
%These functions are the gui controls for the filebrowser function 

feval(varargin{:});
%--------------------------------------------
function CloseRequestFnc(handles)

handles.filename = '';
handles.path = '';
guidata(handles.fighandle, handles);
uiresume(handles.fighandle);
%--------------------------------------------
function cancelbutton_Callback(handles)

handles.filename = '';
handles.path = '';
guidata(handles.fighandle, handles);
uiresume(handles.fighandle);
%----------------------------------------------
function okbutton_Callback(handles)
%the user clicked the ok button.  We do some checks before we allow the
%figure to close
value = get(handles.listbox1,'Value');
if ((handles.dirstruct(value).isdir == 1) & (strcmp(handles.filename,'')))
    cd(handles.dirstruct(value).name);
    set(handles.pathedit,'String',pwd);
    pathedit_Callback(handles);
    return
end
    
if ((handles.filemode==2)&(exist(handles.filename,'file')~=2))
    errordlg(['Filename ',handles.filename,' does not exist.']);
    return
end
if ((handles.filemode==1)&(strcmp(handles.filename,'')))
    errordlg('Please choose a name for the saved file.');
    return
end
if ((handles.filemode==1)&(exist(handles.filename,'file')==2))
    reply = questdlg(['Overwrite ',handles.filename,'?'],'Overwrite?','Cancel','Ok','Ok');
    if (strcmp(reply,'Cancel'))
        return
    end
end

%all checks cleared, so we release the hold on the figure
uiresume(handles.fighandle);
%---------------------------------------------
function listbox1_Callback(handles)
%the user clicked an item in the listbox

value = get(handles.listbox1,'Value');
clicktype = get(handles.fighandle,'SelectionType');
dirlist = handles.dirlist;
dirstruct = handles.dirstruct;
if (strcmp(clicktype,'open') & (handles.dirstruct(value).isdir == 1)) %the user double-clicked on a directory, so we should change the directory
    
    cd(handles.dirstruct(value).name);
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
    diroutput = dir(handles.filter);
    for a = 1:length(diroutput)  
        if ~diroutput(a).isdir     
            dirlist{count} = diroutput(a).name;
            dirstruct(count) = diroutput(a);
            count = count+1;
        end
    end
    handles.filename = '';
    handles.path = pwd;
    set(handles.fileedit,'String','');
    set(handles.listbox1,'Value',1);
    set(handles.listbox1,'String',dirlist);
    set(handles.pathedit,'String',pwd);
end
if (strcmp(clicktype,'normal') & (handles.dirstruct(value).isdir == 0)) %the user single-clicked on a file
    handles.filename = handles.dirstruct(value).name;
    set(handles.fileedit,'String',handles.dirstruct(value).name);
end
if (strcmp(clicktype,'open') & (handles.dirstruct(value).isdir == 0)) %the user double-clicked a file, activate the ok button
    handles.filename = handles.dirstruct(value).name;
    set(handles.fileedit,'String',handles.dirstruct(value).name);
    guidata(handles.fighandle, handles);
    okbutton_Callback(handles)
end
if (strcmp(clicktype,'normal') & (handles.dirstruct(value).isdir == 1))
    handles.filename = '';
    set(handles.fileedit,'String','');
end
if ((handles.dirstruct(value).isdir == -1))
    handles.filename = '';
    set(handles.fileedit,'String','');
end

handles.dirlist = dirlist;
handles.dirstruct = dirstruct;
guidata(handles.fighandle, handles);
%-------------------------------------------------------
function pathedit_Callback(handles)
%the user entered a path manually, so we should change the directory
newpath = get(handles.pathedit,'String');
try
    cd(newpath);
catch
    set(handles.pathedit,'String',pwd);
end
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
diroutput = dir(handles.filter);
for a = 1:length(diroutput)  
    if ~diroutput(a).isdir     
        dirlist{count} = diroutput(a).name;
        dirstruct(count) = diroutput(a);
        count = count+1;
    end
end
set(handles.listbox1,'Value',1);
handles.path = pwd;
handles.dirlist = dirlist;
handles.dirstruct = dirstruct;
set(handles.listbox1,'String',dirlist);
guidata(handles.fighandle, handles);
%-------------------------------------------
function fileedit_Callback(handles)
%the user entered a filename

handles.filename = get(handles.fileedit,'String');
guidata(handles.fighandle, handles);
%----------------------------------------------
function filteredit_Callback(handles)
%the user changed the filter

handles.filter = get(handles.filteredit,'String');
guidata(handles.fighandle, handles);
pathedit_Callback(handles);