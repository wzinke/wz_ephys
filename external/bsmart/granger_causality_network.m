function varargout = granger_causality_network(varargin)
% GRANGER_CAUSALITY_NETWORK M-file for granger_causality_network.fig
%      GRANGER_CAUSALITY_NETWORK, by itself, creates a new GRANGER_CAUSALITY_NETWORK or raises the existing
%      singleton*.
%
%      H = GRANGER_CAUSALITY_NETWORK returns the handle to a new GRANGER_CAUSALITY_NETWORK or the handle to
%      the existing singleton*.
%
%      GRANGER_CAUSALITY_NETWORK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GRANGER_CAUSALITY_NETWORK.M with the given input arguments.
%
%      GRANGER_CAUSALITY_NETWORK('Property','Value',...) creates a new GRANGER_CAUSALITY_NETWORK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before granger_causality_network_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to granger_causality_network_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.3$ $Date: 15-Sep-2007 19:21:43$
% SHIS UT-Houston, Houston, TX 77030, USA.
%
% Lei Xu, Hualou Liang

% Edit the above text to modify the response to help granger_causality_network

% Last Modified by GUIDE v2.5 15-Sep-2007 21:20:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @granger_causality_network_OpeningFcn, ...
                   'gui_OutputFcn',  @granger_causality_network_OutputFcn, ...
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


% --- Executes just before granger_causality_network is made visible.
function granger_causality_network_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to granger_causality_network (see VARARGIN)

% Choose default command line output for granger_causality_network
handles.output = hObject;

% Update handles structure
handles.fxyfilepath = '';
handles.fyxfilepath = '';
handles.locfilepath = '';
guidata(hObject, handles);

% UIWAIT makes granger_causality_network wait for user response (see UIRESUME)
% uiwait(handles.gc_network);


% --- Outputs from this function are returned to the command line.
function varargout = granger_causality_network_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in chose_Fxy.
function chose_Fxy_Callback(hObject, eventdata, handles)
% hObject    handle to chose_Fxy (see GCBO)
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
else
    handles.fxyfilepath = fullfile(pathname,filename);
end
guidata(hObject,handles);
set(handles.fFxy_show,'String',filename);


function fFxy_show_Callback(hObject, eventdata, handles)
% hObject    handle to fFxy_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fFxy_show as text
%        str2double(get(hObject,'String')) returns contents of fFxy_show as a double


% --- Executes during object creation, after setting all properties.
function fFxy_show_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fFxy_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function fstart_Callback(hObject, eventdata, handles)
% hObject    handle to fstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fstart as text
%        str2double(get(hObject,'String')) returns contents of fstart as a double


% --- Executes during object creation, after setting all properties.
function fstart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fend_Callback(hObject, eventdata, handles)
% hObject    handle to fend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fend as text
%        str2double(get(hObject,'String')) returns contents of fend as a double


% --- Executes during object creation, after setting all properties.
function fend_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function win_Callback(hObject, eventdata, handles)
% hObject    handle to win (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of win as text
%        str2double(get(hObject,'String')) returns contents of win as a double


% --- Executes during object creation, after setting all properties.
function win_CreateFcn(hObject, eventdata, handles)
% hObject    handle to win (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

                delete(handles.gc_network);

% --- Executes on button press in show_net.
function show_net_Callback(hObject, eventdata, handles)
% hObject    handle to show_net (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Fxy file
if isequal(handles.fxyfilepath,'')
    errordlg('File not found','File Error');
    return
end
a    = load(handles.fxyfilepath);
name = fieldnames(a);
d    = name{1};
fxy  = getfield(a,d);

% Fyx file
if isequal(handles.fyxfilepath,'')
    errordlg('File not found','File Error');
    return
end
a    = load(handles.fyxfilepath);
name = fieldnames(a);
d    = name{1};
fyx  = getfield(a,d);

% location file
if isequal(handles.locfilepath,'')
    errordlg('File not found','File Error');
    return
end
a    = load(handles.locfilepath);
name = fieldnames(a);
d    = name{1};
loc  = getfield(a,d);

fre1_string=get(handles.fstart,'String');
if fre1_string==' '
    errordlg('please input correct frequency','parameter lost');
    return
end
fstart = str2double(fre1_string);
if (fstart <= 0)
    errordlg('please input correct frequency','parameter lost');
    return
end

fre2_string=get(handles.fend,'String');
if fre2_string==' '
    errordlg('please input correct frequency','parameter lost');
    return
end
fend = str2double(fre2_string);
if (fend <= 0)
    errordlg('please input correct frequency','parameter lost');
    return
end
if (fstart > fend)
    errordlg('please input correct start and end frequency','parameter lost');
    return
end

time_string = get(handles.win,'String');
if time_string==' '
       errordlg('please input window','parameter lost');
       return
end
win = str2double(time_string);
if win < 0
    errordlg('please input correct window','parameter lost');
    return
end

thre_string = get(handles.thres,'String');
if thre_string==' '
    errordlg('please input threshold','parameter lost');
    return
end
thre = str2double(thre_string);
if thre < 0
    errordlg('please input correct threshold','parameter lost');
    return
end

% set chan of interest 
s1  = size(fxy);
c   = s1(1);
nch = (1+sqrt(1+8*c))/2;    % number of channels
chint_string = get(handles.ch_int,'String');
if isempty(chint_string)
    chan = 1:nch;
else
    chan = str2num(chint_string); %#ok<ST2NM>
end

% core: shaw gc network
ganetwork(fxy,fyx,loc,thre,win,fstart,fend,chan)


function thres_Callback(hObject, eventdata, handles)
% hObject    handle to thres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thres as text
%        str2double(get(hObject,'String')) returns contents of thres as a double


% --- Executes during object creation, after setting all properties.
function thres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in choose_Fyx.
function choose_Fyx_Callback(hObject, eventdata, handles)
% hObject    handle to choose_Fyx (see GCBO)
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
else
    handles.fyxfilepath = fullfile(pathname,filename);
end
guidata(hObject,handles);
set(handles.fFyx_show,'String',filename);

function fFyx_show_Callback(hObject, eventdata, handles)
% hObject    handle to fFyx_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fFyx_show as text
%        str2double(get(hObject,'String')) returns contents of fFyx_show as a double


% --- Executes during object creation, after setting all properties.
function fFyx_show_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fFyx_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in choose_loc.
function choose_loc_Callback(hObject, eventdata, handles)
% hObject    handle to choose_loc (see GCBO)
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
else
    handles.locfilepath = fullfile(pathname,filename);
end
guidata(hObject,handles);
set(handles.floc_show,'String',filename);


function floc_show_Callback(hObject, eventdata, handles)
% hObject    handle to floc_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of floc_show as text
%        str2double(get(hObject,'String')) returns contents of floc_show as a double


% --- Executes during object creation, after setting all properties.
function floc_show_CreateFcn(hObject, eventdata, handles)
% hObject    handle to floc_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ch_int_Callback(hObject, eventdata, handles)
% hObject    handle to ch_int (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch_int as text
%        str2double(get(hObject,'String')) returns contents of ch_int as a double


% --- Executes during object creation, after setting all properties.
function ch_int_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch_int (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


