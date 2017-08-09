function varargout = preprocessing(varargin)
% PREPROCESSING M-file for preprocessing.fig
%      PREPROCESSING, by itself, creates a new PREPROCESSING or raises the existing
%      singleton*.
%
%      H = PREPROCESSING returns the handle to a new PREPROCESSING or the handle to
%      the existing singleton*.
%
%      PREPROCESSING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PREPROCESSING.M with the given input arguments.
%
%      PREPROCESSING('Property','Value',...) creates a new PREPROCESSING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before preprocessing_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to preprocessing_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help preprocessing

% Last Modified by GUIDE v2.5 14-Sep-2007 15:54:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @preprocessing_OpeningFcn, ...
                   'gui_OutputFcn',  @preprocessing_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before preprocessing is made visible.
function preprocessing_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to preprocessing (see VARARGIN)

% Choose default command line output for preprocessing
handles.output = hObject;

% Update handles structure
handles.nchannels=0;
handles.ntrails=0;
handles.npoints=0;
handles.filepath='';
guidata(hObject, handles);

% UIWAIT makes preprocessing wait for user response (see UIRESUME)
% uiwait(handles.preprocessing);


% --- Outputs from this function are returned to the command line.
function varargout = preprocessing_OutputFcn(hObject, eventdata, handles)  %#ok<INUSL>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in choose_mat_file.
function choose_mat_file_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to choose_mat_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile( ...
    { '*.mat','MAT-files (*.mat)'; ...
    '*.*',  'All Files (*.*)'},...
    'Pick a file');

if isequal(filename,0) || isequal(pathname,0)
    filename = 'No file selected';
else
    handles.filepath=fullfile(pathname,filename);
end

set(handles.fname_show,'String',filename);

guidata(hObject,handles);


function fname_show_Callback(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
% hObject    handle to fname_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fname_show as text
%        str2double(get(hObject,'String')) returns contents of fname_show as a double


% --- Executes during object creation, after setting all properties.
function fname_show_CreateFcn(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
% hObject    handle to fname_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.


% --- Executes on button press in sube.
function sube_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to sube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(handles.filepath,'')
    errordlg('File not found','File Error');
    return;
end

oldfullname = handles.filepath;
data  = load(oldfullname);
name  = fieldnames(data);
d     = name{1};
data  = getfield(data,d); %#ok<GFLD>

% core: substruct ensemble mean
dat = pre_sube(data); %#ok<NASGU>

% save file
[pathstr, suggestname] = fileparts(oldfullname);
[filename,pathname] = uiputfile( ...
    {'*.mat','MAT-files (*.mat)';...
    '*.*',  'All Files (*.*)'},...
    'Save as',...
    suggestname);

if isequal(filename,0) || isequal(pathname,0)
else
    filepath = fullfile(pathname,filename);
    save(filepath,'dat');
end%if


% --- Executes on button press in sube_divs.
function sube_divs_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to sube_divs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(handles.filepath,'')
    errordlg('File not found','File Error');
    return;
end

oldfullname = handles.filepath;
data    = load(oldfullname);
name    = fieldnames(data);
d       = name{1};
data    = getfield(data,d); %#ok<GFLD>

% core: sube_divs
dat = pre_sube_divs(data); %#ok<NASGU>

% save file
[pathstr, suggestname] = fileparts(oldfullname);
[filename,pathname] = uiputfile( ...
    {'*.mat','MAT-files (*.mat)';...
    '*.*',  'All Files (*.*)'},...
    'Save as',...
    suggestname);

if isequal(filename,0) || isequal(pathname,0)
else
    filepath = fullfile(pathname,filename);
    save(filepath,'dat');
end%if


% --- Executes on button press in subt.
function subt_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to subt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(handles.filepath,'')
    errordlg('File not found','File Error');
    return
end

oldfullname = handles.filepath;
data    = load(oldfullname);
name    = fieldnames(data);
d       = name{1};
data    = getfield(data,d); %#ok<GFLD>

% core: pre_subt
dat = pre_subt(data); %#ok<NASGU>

% save file
[pathstr, suggestname] = fileparts(oldfullname);
[filename,pathname] = uiputfile( ...
    {'*.mat','MAT-files (*.mat)';...
    '*.*',  'All Files (*.*)'},...
    'Save as',...
    suggestname);

if isequal(filename,0) || isequal(pathname,0)
else
    filepath = fullfile(pathname,filename);
    save(filepath,'dat');
end%if


% --- Executes on button press in subt_divs.
function subt_divs_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to subt_divs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(handles.filepath,'')
    errordlg('File not found','File Error');
    return
end

oldfullname = handles.filepath;
data    = load(oldfullname);
name    = fieldnames(data);
d       = name{1};
data    = getfield(data,d); %#ok<GFLD>

% core: pre_subt_divs
dat = pre_subt_divs(data); %#ok<NASGU>

% save file
[pathstr, suggestname] = fileparts(oldfullname);
[filename,pathname] = uiputfile( ...
    {'*.mat','MAT-files (*.mat)';...
    '*.*',  'All Files (*.*)'},...
    'Save as',...
    suggestname);

if isequal(filename,0) || isequal(pathname,0)
else
    filepath = fullfile(pathname,filename);
    save(filepath,'dat');
end%if


% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.preprocessing);

