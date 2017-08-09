function varargout = writedata(varargin)
% WRITEDATA M-file for writedata.fig
%      WRITEDATA, by itself, creates a new WRITEDATA or raises the existing
%      singleton*.
%
%      H = WRITEDATA returns the handle to a new WRITEDATA or the handle to
%      the existing singleton*.
%
%      WRITEDATA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WRITEDATA.M with the given input arguments.
%
%      WRITEDATA('Property','Value',...) creates a new WRITEDATA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before writedata_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to writedata_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help writedata

% Last Modified by GUIDE v2.5 14-Sep-2007 14:45:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @writedata_OpeningFcn, ...
                   'gui_OutputFcn',  @writedata_OutputFcn, ...
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


% --- Executes just before writedata is made visible.
function writedata_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to writedata (see VARARGIN)

% Choose default command line output for writedata
handles.output = hObject;

% Update handles structure
handles.filepath='';
guidata(hObject, handles);

% UIWAIT makes writedata wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = writedata_OutputFcn(hObject, eventdata, handles)  %#ok<INUSL>
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

% select a file
[filename, pathname] = uigetfile( ...
    {'*.mat','MAT-files (*.mat)'; ...
    '*.*',  'All Files (*.*)'},...
    'Pick a file');

if isequal(filename,0) || isequal(pathname,0)
    filename = 'No file selected';
else
    handles.filepath = fullfile(pathname,filename);
end

% show file name
set(handles.fname_show,'String',filename);

% save data structure
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



% --- Executes on button press in export_data.
function export_data_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to export_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(handles.filepath,'')
    errordlg('File not found','File Error');
    return
end

[file,path] = uiputfile('*.bin','Save as');

if isequal(file,0) || isequal(path,0)
else
    oldfile = handles.filepath;
    newfile = fullfile(path,file);
    c = load(oldfile);
    name = fieldnames(c);
    d = name{1};
    dat = getfield(c,d); %#ok<GFLD>
    writedat(newfile,dat);
end

% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.export);

