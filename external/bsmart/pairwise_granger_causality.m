function varargout = pairwise_granger_causality(varargin)
% PAIRWISE_GRANGER_CAUSALITY M-file for pairwise_granger_causality.fig
%      PAIRWISE_GRANGER_CAUSALITY, by itself, creates a new PAIRWISE_GRANGER_CAUSALITY or raises the existing
%      singleton*.
%
%      H = PAIRWISE_GRANGER_CAUSALITY returns the handle to a new PAIRWISE_GRANGER_CAUSALITY or the handle to
%      the existing singleton*.
%
%      PAIRWISE_GRANGER_CAUSALITY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PAIRWISE_GRANGER_CAUSALITY.M with the given input arguments.
%
%      PAIRWISE_GRANGER_CAUSALITY('Property','Value',...) creates a new PAIRWISE_GRANGER_CAUSALITY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pairwise_granger_causality_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pairwise_granger_causality_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help pairwise_granger_causality

% Last Modified by GUIDE v2.5 15-Sep-2007 17:11:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pairwise_granger_causality_OpeningFcn, ...
                   'gui_OutputFcn',  @pairwise_granger_causality_OutputFcn, ...
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


% --- Executes just before pairwise_granger_causality is made visible.
function pairwise_granger_causality_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pairwise_granger_causality (see VARARGIN)

% Choose default command line output for pairwise_granger_causality
handles.output = hObject;

% Update handles structure
handles.filepath = '';
guidata(hObject, handles);

% UIWAIT makes pairwise_granger_causality wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pairwise_granger_causality_OutputFcn(hObject, eventdata, handles) 
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
    '*.m;*.fig;*.mat;*.mdl','MATLAB Files (*.m,*.fig,*.mat,*.mdl)';...
    '*.m',  'M-files (*.m)'; ...
    '*.fig','Figures (*.fig)'; ...
    '*.mdl','Models (*.mdl)'; ...
    '*.*',  'All Files (*.*)'; ...
    '*.dm6','steven file'},...
    'Pick a file');
if isequal(filename,0) || isequal(pathname,0)
    filepath = handles.filepath;
    if isempty(filepath)
        filename = 'No file selected...';
    else
        [path,name,ext] = fileparts(filepath);
        filename = [name,ext];
    end%of
    %disp('User selected Cancel')
else
    handles.filepath = fullfile(pathname,filename);
end
guidata(hObject,handles);
set(handles.fname_show,'String',filename);


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
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
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
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function fs_Callback(hObject, eventdata, handles)
% hObject    handle to fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fs as text
%        str2double(get(hObject,'String')) returns contents of fs as a double


% --- Executes during object creation, after setting all properties.
function fs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function frang_Callback(hObject, eventdata, handles)
% hObject    handle to frang (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frang as text
%        str2double(get(hObject,'String')) returns contents of frang as a double


% --- Executes during object creation, after setting all properties.
function frang_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frang (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in gc_mov.
function gc_mov_Callback(hObject, eventdata, handles)
% hObject    handle to gc_mov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% parameter input
if isequal(handles.filepath,'')
    errordlg('File not found','File Error');
    return
end
data  = load(handles.filepath);
name  = fieldnames(data);
d     = name{1};
dat   = getfield(data,d);

order_string  = get(handles.morder,'String');
if order_string==' '
    errordlg('please input the model order','parameter lost');
    return
end
morder = str2double(order_string);

fs_string     = get(handles.fs,'String');
if fs_string==' '
    errordlg('please input sample frequency','parameter lost');
    return
end
fs = str2double(fs_string);

start_string  = get(handles.spts,'String');
if start_string==' '
    errordlg('please input the starting position','parameter lost');
    return
end
spts  = str2double(start_string);

window_string = get(handles.winlen,'String');
if window_string==' '
    errordlg('please input the window','parameter lost');
    return
end
winlen = str2double(window_string);

freq_string   = get(handles.frang,'String');
if freq_string==' '
    errordlg('please input a vector of frequencies of interest, usually freq=0:fs/2','parameter lost');
    return
end
frang = str2num(freq_string);  %#ok<ST2NM>

epts_string = get(handles.epts,'String');
if freq_string==' '
    errordlg('please input ending point','parameter lost');
    return
end
epts = str2double(epts_string);

% core:
[Fxy,Fyx] = mov_bi_ga(dat,spts,epts,winlen,morder,fs,frang); %#ok<NASGU>

% save
suggestname = 'Fxy';
[filename,pathname] = uiputfile( ...
    {'*.mat','MAT-files (*.mat)';...
    '*.*',  'All Files (*.*)'},...
    'Save as',...
    suggestname);
if isequal(filename,0) || isequal(pathname,0)
else
    filepath = fullfile(pathname,filename);
    save(filepath,'Fxy');
end%if

suggestname = 'Fyx';
[filename,pathname] = uiputfile( ...
    {'*.mat','MAT-files (*.mat)';...
    '*.*',  'All Files (*.*)'},...
    'Save as',...
    suggestname);
if isequal(filename,0) || isequal(pathname,0)
else
    filepath = fullfile(pathname,filename);
    save(filepath,'Fyx');
end%if


% --- Executes on button press in gc_fix.
function gc_fix_Callback(hObject, eventdata, handles)
% hObject    handle to gc_fix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% parameter input
if isequal(handles.filepath,'')
    errordlg('File not found','File Error');
    return
end
data  = load(handles.filepath);
name  = fieldnames(data);
d     = name{1};
dat   = getfield(data,d);

order_string  = get(handles.morder,'String');
if order_string==' '
    errordlg('please input the model order','parameter lost');
    return
end
morder = str2double(order_string);

fs_string     = get(handles.fs,'String');
if fs_string==' '
    errordlg('please input sample frequency','parameter lost');
    return
end
fs = str2double(fs_string);

start_string  = get(handles.spts,'String');
if start_string==' '
    errordlg('please input the starting position','parameter lost');
    return
end
spts  = str2double(start_string);

window_string = get(handles.winlen,'String');
if window_string==' '
    errordlg('please input the window','parameter lost');
    return
end
winlen = str2double(window_string);

freq_string   = get(handles.frang,'String');
if freq_string==' '
    errordlg('please input a vector of frequencies of interest, usually freq=0:fs/2','parameter lost');
    return
end
frang = str2num(freq_string); %#ok<ST2NM>

% core:
h = msgbox('Please wait...','Granger causality','help');
[Fx2y,Fy2x]= one_bi_ga(dat,spts,winlen,morder,fs,frang); %#ok<NASGU>
close(h);

% save
suggestname = 'Fx2y';
[filename,pathname] = uiputfile( ...
    {'*.mat','MAT-files (*.mat)';...
    '*.*',  'All Files (*.*)'},...
    'Save as',...
    suggestname);
if isequal(filename,0) || isequal(pathname,0)
else
    filepath = fullfile(pathname,filename);
    save(filepath,'Fx2y');
end%if

suggestname = 'Fy2x';
[filename,pathname] = uiputfile( ...
    {'*.mat','MAT-files (*.mat)';...
    '*.*',  'All Files (*.*)'},...
    'Save as',...
    suggestname);
if isequal(filename,0) || isequal(pathname,0)
else
    filepath = fullfile(pathname,filename);
    save(filepath,'Fy2x');
end%if


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
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);


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


