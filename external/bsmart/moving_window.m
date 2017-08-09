function varargout = moving_window(varargin)
% MOVING_WINDOW M-file for moving_window.fig
%      MOVING_WINDOW, by itself, creates a new MOVING_WINDOW or raises the existing
%      singleton*.
%
%      H = MOVING_WINDOW returns the handle to a new MOVING_WINDOW or the handle to
%      the existing singleton*.
%
%      MOVING_WINDOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MOVING_WINDOW.M with the given input arguments.
%
%      MOVING_WINDOW('Property','Value',...) creates a new MOVING_WINDOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before moving_window_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to moving_window_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

%   Copyright 2006-2007 BSMART group
%   $Revision: 0.2 $  $Date: 03-Jul-2007 14:34:59 $
%   SHIS UT-Houston, Houston, TX 77030, USA.
%   Coded by: X. Lei, R.J. Cui

% Edit the above text to modify the response to help moving_window

% Last Modified by GUIDE v2.5 14-Sep-2007 21:55:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @moving_window_OpeningFcn, ...
                   'gui_OutputFcn',  @moving_window_OutputFcn, ...
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


% --- Executes just before moving_window is made visible.
function moving_window_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to moving_window (see VARARGIN)

% Choose default command line output for moving_window
handles.output = hObject;

% Update handles structure
handles.nchannels=0;
handles.ntrails=0;
handles.npoints=0;
handles.window=0;
handles.order=0;
handles.filepath='';
guidata(hObject, handles);

% UIWAIT makes moving_window wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = moving_window_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in choose_mat_file.
function choose_mat_file_Callback(hObject, eventdata, handles)
% hObject    handle to choose_mat_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
{'*.mat','MAT-files (*.mat)'; ...
'*.m;*.fig;*.mat;*.mdl','MATLAB Files (*.m,*.fig,*.mat,*.mdl)';
   '*.m',  'M-files (*.m)'; ...
   '*.fig','Figures (*.fig)'; ...
   '*.mdl','Models (*.mdl)'; ...
   '*.*',  'All Files (*.*)'; ...
   '*.dm6','steven file'},...
   'Pick a file');
if isequal(filename,0) | isequal(pathname,0)
   %disp('User selected Cancel')
else
handles.filepath=fullfile(pathname,filename);
end
guidata(hObject,handles);
set(handles.fname_show,'String',handles.filepath);


function fname_show_Callback(hObject, eventdata, handles)
% hObject    handle to fname_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fname_show as text
%        str2double(get(hObject,'String')) returns contents of fname_show as a double


% --- Executes during object creation, after setting all properties.
function fname_show_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fname_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in create_model.
function create_model_Callback(hObject, eventdata, handles)  %#ok<DEFNU,INUSL>
% hObject    handle to create_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% process inputs
start_string = get(handles.spts,'String');
end_string   = get(handles.epts,'String');
window_string= get(handles.winlen,'String');
order_string = get(handles.morder,'String');
if isequal(handles.filepath,'')
    errordlg('File not found','File Error');
    return;
end

data = load(handles.filepath);
name = fieldnames(data);
d    = name{1,1};        % only get the 1st data struct
dat  = getfield(data,d); %#ok<GFLD>
order  = str2double(order_string);
spts   = str2double(start_string);
epts   = str2double(end_string);
winlen = str2double(window_string);

% core: moving window multivariate model
h = msgbox('Please waiting...','Moving window multivariate model','help');
[A,Ve] = mov_mul_model(dat,order,spts,epts,winlen); %#ok<NASGU>
close(h);

% save results
% (1) AR coefficients
suggestname = 'ARC_mul';
[filename,pathname] = uiputfile( ...
    {'*.mat','MAT-files (*.mat)';...
    '*.*',  'All Files (*.*)'},...
    'Save AR Coefficients',...
    suggestname);

if isequal(filename,0) || isequal(pathname,0)
else
    filepath = fullfile(pathname,filename);
    save(filepath,'A');
end%if
% (2) Noise covariance matrix
suggestname = 'NCM_mul';
[filename,pathname] = uiputfile( ...
    {'*.mat','MAT-files (*.mat)';...
    '*.*',  'All Files (*.*)'},...
    'Save Covariance Matrix of Noise Vector',...
    suggestname);

if isequal(filename,0) || isequal(pathname,0)
else
    filepath = fullfile(pathname,filename);
    save(filepath,'Ve');
end%if


% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
                delete(handles.figure1);


function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function winlen_Callback(hObject, eventdata, handles)
% hObject    handle to winlen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of winlen as text
%        str2double(get(hObject,'String')) returns contents of winlen as a double


% --- Executes during object creation, after setting all properties.
function winlen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to winlen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function morder_Callback(hObject, eventdata, handles)
% hObject    handle to morder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of morder as text
%        str2double(get(hObject,'String')) returns contents of morder as a double


% --- Executes during object creation, after setting all properties.
function morder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to morder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function epts_Callback(hObject, eventdata, handles)
% hObject    handle to epts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of epts as text
%        str2double(get(hObject,'String')) returns contents of epts as a double


% --- Executes during object creation, after setting all properties.
function epts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function spts_Callback(hObject, eventdata, handles)
% hObject    handle to spts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spts as text
%        str2double(get(hObject,'String')) returns contents of spts as a double


% --- Executes during object creation, after setting all properties.
function spts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% [EOF]