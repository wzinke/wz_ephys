function varargout = power_view(varargin)
% POWER_VIEW M-file for power_view.fig
%      POWER_VIEW, by itself, creates a new POWER_VIEW or raises the existing
%      singleton*.
%
%      H = POWER_VIEW returns the handle to a new POWER_VIEW or the handle to
%      the existing singleton*.
%
%      POWER_VIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POWER_VIEW.M with the given input arguments.
%
%      POWER_VIEW('Property','Value',...) creates a new POWER_VIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before power_view_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to power_view_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright (c) 2006-2007 BSMART Group
% by Richard Cui
% $Revision: 0.3$ $Date: 19-Sep-2007 16:07:55$
% SHIS UT-Houston, Houston, TX 77030, USA.
% 
% Lei Xu, Hualou Liang

% Edit the above text to modify the response to help power_view

% Last Modified by GUIDE v2.5 15-Sep-2007 23:29:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @power_view_OpeningFcn, ...
                   'gui_OutputFcn',  @power_view_OutputFcn, ...
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


% --- Executes just before power_view is made visible.
function power_view_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to power_view (see VARARGIN)

% Choose default command line output for power_view
handles.output = hObject;

% Update handles structure
handles.filepath='';
handles.channel=0;
guidata(hObject, handles);

% UIWAIT makes power_view wait for user response (see UIRESUME)
% uiwait(handles.power_view);


% --- Outputs from this function are returned to the command line.
function varargout = power_view_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function chan_Callback(hObject, eventdata, handles)
% hObject    handle to chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of chan as text
%        str2double(get(hObject,'String')) returns contents of chan as a double


% --- Executes during object creation, after setting all properties.
function chan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in plot_timefreq.
function plot_timefreq_Callback(hObject, eventdata, handles)
% hObject    handle to plot_timefreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% power
if isequal(handles.filepath,'')
    errordlg('File not found','File Error');
    return
end

data   = load(handles.filepath);
name   = fieldnames(data);
d      = name{1};
power  = getfield(data,d);

% parameters
nwin   = size(power,1);     % number of windows
ncha   = size(power,3);     % number of channels

% sampling rate
fs_string  = get(handles.fs,'String');
if fs_string==' '
        errordlg('please input the sampling rate','parameter lost');
        return
end
fs = str2double(fs_string);
if fs <= 0
    errordlg('Incorrect sampling rate','Parameter Error');
    return
end%if

% specified channel
chan_string = get(handles.chan,'String');
if chan_string==' '
       errordlg('please input channel number','parameter lost');
       return
end
chan = str2double(chan_string);
if chan <= 0 || chan > ncha
    errordlg('Invalid channel number','Parameter Error');
    return
end%if

% core: power view
po_view(power,fs,chan);


% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 delete(handles.power_view);


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
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in plot_spect.
function plot_spect_Callback(hObject, eventdata, handles)
% hObject    handle to plot_spect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% power
if isequal(handles.filepath,'')
    errordlg('File not found','File Error');
    return
end

data   = load(handles.filepath);
name   = fieldnames(data);
d      = name{1};
power  = getfield(data,d);
% parameters
nwin   = size(power,1);     % number of windows
ncha   = size(power,3);     % number of channels

% sampling rate
fs_string  = get(handles.fs,'String');
if fs_string==' '
        errordlg('please input the sampling rate','parameter lost');
        return
end
fs = str2double(fs_string);
if fs <= 0
    errordlg('Incorrect sampling rate','Parameter Error');
    return
end%if

% specified channel
chan_string = get(handles.chan,'String');
if chan_string==' '
       errordlg('please input channel number','parameter lost');
       return
end
chan = str2double(chan_string);
if chan <= 0 || chan > ncha
    errordlg('Invalid channel number','Parameter Error');
    return
end%if

% specified window
time_string  = get(handles.win,'String');
if time_string==' '
    errordlg('please input window','parameter lost');
    return
end
win = str2double(time_string);
if win <= 0 || win > nwin
    errordlg('Please input correct window','Parameter Error');
    return
end%if

% core: power view
po_view(power,fs,chan,win);

   
% --- Executes on button press in choose_coh.
function choose_coh_Callback(hObject, eventdata, handles)
% hObject    handle to choose_coh (see GCBO)
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
if isequal(filename,0) | isequal(pathname,0)
else
    handles.filepath=fullfile(pathname,filename);
end
guidata(hObject,handles);
set(handles.fcoh_show,'String',handles.filepath);


function fcoh_show_Callback(hObject, eventdata, handles)
% hObject    handle to fcoh_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fcoh_show as text
%        str2double(get(hObject,'String')) returns contents of fcoh_show as a double


% --- Executes during object creation, after setting all properties.
function fcoh_show_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fcoh_show (see GCBO)
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


