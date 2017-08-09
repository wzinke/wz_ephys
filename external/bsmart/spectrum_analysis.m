function varargout = spectrum_analysis(varargin)
% SPECTRUM_ANALYSIS M-file for spectrum_analysis.fig
%      SPECTRUM_ANALYSIS, by itself, creates a new SPECTRUM_ANALYSIS or raises the existing
%      singleton*.
%
%      H = SPECTRUM_ANALYSIS returns the handle to a new SPECTRUM_ANALYSIS or the handle to
%      the existing singleton*.
%
%      SPECTRUM_ANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPECTRUM_ANALYSIS.M with the given input arguments.
%
%      SPECTRUM_ANALYSIS('Property','Value',...) creates a new SPECTRUM_ANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before spectrum_analysis_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to spectrum_analysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright (c) 2006-2007 BSMART Group
% by Richard Cui
% $Revision: 0.3$ $Date: 19-Sep-2007 16:18:17$
% SHIS UT-Houston, Houston, TX 77030, USA.
% 
% Lei Xu, Hualou Liang

% Edit the above text to modify the response to help spectrum_analysis

% Last Modified by GUIDE v2.5 15-Sep-2007 13:31:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @spectrum_analysis_OpeningFcn, ...
                   'gui_OutputFcn',  @spectrum_analysis_OutputFcn, ...
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


% --- Executes just before spectrum_analysis is made visible.
function spectrum_analysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to spectrum_analysis (see VARARGIN)

% Choose default command line output for spectrum_analysis
handles.output = hObject;

% Update handles structure
handles.arcfilepath = '';
handles.arnfilepath = '';
guidata(hObject, handles);

% UIWAIT makes spectrum_analysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = spectrum_analysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function choose_arc_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile( ...
    {'*.*',  'All Files (*.*)'; ...
    '*.m;*.fig;*.mat;*.mdl','MATLAB Files (*.m,*.fig,*.mat,*.mdl)';...
    '*.m',  'M-files (*.m)'; ...
    '*.fig','Figures (*.fig)'; ...
    '*.mat','MAT-files (*.mat)'; ...
    '*.mdl','Models (*.mdl)'; ...
    '*.dm6','steven file'},...
    'Pick a file');
if isequal(filename,0) || isequal(pathname,0)
else
    handles.arcfilepath=fullfile(pathname,filename);
end
guidata(hObject,handles);
set(handles.farc_show,'String',handles.arcfilepath);

% --- Executes on button press in choose_ncm.
function choose_ncm_Callback(hObject, eventdata, handles)
% hObject    handle to choose_ncm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {'*.*',  'All Files (*.*)'; ...
    '*.m;*.fig;*.mat;*.mdl','MATLAB Files (*.m,*.fig,*.mat,*.mdl)';...
    '*.m',  'M-files (*.m)'; ...
    '*.fig','Figures (*.fig)'; ...
    '*.mat','MAT-files (*.mat)'; ...
    '*.mdl','Models (*.mdl)'; ...
    '*.dm6','steven file'},...
    'Pick a file');
if isequal(filename,0) || isequal(pathname,0)
else
    handles.arnfilepath=fullfile(pathname,filename);
end
guidata(hObject,handles);
set(handles.fncm_show,'String',handles.arnfilepath);



function farc_show_Callback(hObject, eventdata, handles)
% hObject    handle to farc_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of farc_show as text
%        str2double(get(hObject,'String')) returns contents of farc_show as a double


% --- Executes during object creation, after setting all properties.
function farc_show_CreateFcn(hObject, eventdata, handles)
% hObject    handle to farc_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fncm_show_Callback(hObject, eventdata, handles)
% hObject    handle to fncm_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fncm_show as text
%        str2double(get(hObject,'String')) returns contents of fncm_show as a double


% --- Executes during object creation, after setting all properties.
function fncm_show_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fncm_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
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
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in power.
function power_Callback(hObject, eventdata, handles)
% hObject    handle to power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% check filepath
if isequal(handles.arcfilepath,'')
    errordlg('Coffecient File not found','File Error');
    return;
end
if isequal(handles.arnfilepath,'')
    errordlg('Noise File not found','File Error');
    return;
end

% check parameters
fs_string  = get(handles.fs,'String');
n_string   = get(handles.nbin,'String');
if fs_string==' '
    errordlg('please input sample frequency','parameter lost');
    return
end
if n_string==' '
    errordlg('please input the length of returned spectrum','parameter lost');
    return
end


nbin = str2double(n_string);
fs   = str2double(fs_string);
% arc   
a    = load(handles.arcfilepath);
name = fieldnames(a);
d    = name{1};
arc  = getfield(a,d); %#ok<GFLD>
% make it compatable with Matlab R14
if size(arc,2) == 1
    arc = arc.';    % if arc is a column, change it to row
end%if

% noise covariance matrix
a    = load(handles.arnfilepath);
name = fieldnames(a);
d    = name{1};
ncm  = getfield(a,d); %#ok<GFLD>
if size(ncm,2) == 1
    ncm = ncm.';    % if ncm is a column, change it to row
end%if

% core: autopower multivariate model
power = autopower(arc,ncm,nbin,fs); %#ok<NASGU>

% save
suggestname = 'power_mul';
[filename,pathname] = uiputfile( ...
    {'*.mat','MAT-files (*.mat)';...
    '*.*',  'All Files (*.*)'},...
    'Save as',...
    suggestname);
if isequal(filename,0) || isequal(pathname,0)
else
    filepath = fullfile(pathname,filename);
    save(filepath,'power');
end%if


% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.figure1);


function nbin_Callback(hObject, eventdata, handles)
% hObject    handle to nbin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nbin as text
%        str2double(get(hObject,'String')) returns contents of nbin as a double


% --- Executes during object creation, after setting all properties.
function nbin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nbin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


