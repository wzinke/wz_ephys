function varargout = coherence_view(varargin)
% COHERENCE_VIEW M-file for coherence_view.fig
%      COHERENCE_VIEW, by itself, creates a new COHERENCE_VIEW or raises the existing
%      singleton*.
%
%      H = COHERENCE_VIEW returns the handle to a new COHERENCE_VIEW or the handle to
%      the existing singleton*.
%
%      COHERENCE_VIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COHERENCE_VIEW.M with the given input arguments.
%
%      COHERENCE_VIEW('Property','Value',...) creates a new COHERENCE_VIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before coherence_view_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to coherence_view_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.4$ $Date: 16-Sep-2007 12:10:03$
% SHIS UT-Houston, Houston, TX 77030, USA.
%
% Lei Xu, Hualou Liang


% Edit the above text to modify the response to help coherence_view

% Last Modified by GUIDE v2.5 16-Sep-2007 12:09:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @coherence_view_OpeningFcn, ...
                   'gui_OutputFcn',  @coherence_view_OutputFcn, ...
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


% --- Executes just before coherence_view is made visible.
function coherence_view_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to coherence_view (see VARARGIN)

% Choose default command line output for coherence_view
handles.output = hObject;

% Update handles structure
handles.filepath = '';
guidata(hObject, handles);

% UIWAIT makes coherence_view wait for user response (see UIRESUME)
% uiwait(handles.coh_view);


% --- Outputs from this function are returned to the command line.
function varargout = coherence_view_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function chx_Callback(hObject, eventdata, handles)
% hObject    handle to chx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of chx as text
%        str2double(get(hObject,'String')) returns contents of chx as a double


% --- Executes during object creation, after setting all properties.
function chx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function chy_Callback(hObject, eventdata, handles)
% hObject    handle to chy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of chy as text
%        str2double(get(hObject,'String')) returns contents of chy as a double


% --- Executes during object creation, after setting all properties.
function chy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in plot_tf.
function plot_tf_Callback(hObject, eventdata, handles)
% hObject    handle to plot_tf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% coherence data
if isequal(handles.filepath,'')
    errordlg('File not found','File Error');
    return
end
data   = load(handles.filepath);
name   = fieldnames(data);
d      = name{1};
coh    = getfield(data,d);
% channel number
c = size(coh,3);
nch = (1+sqrt(1+8*c))/2;

% sampling rate
fs_string = get(handles.fs,'String');
if fs_string == ' '
    errordlg('please input sampling rate','parameter lost');
    return
end
fs = str2double(fs_string);
if (fs <= 0)
    errordlg('please input correct sampling rate','parameter lost');
    return
end

% channel X
chx_string = get(handles.chx,'String');
if chx_string == ' '
    errordlg('please input channel number','parameter lost');
    return
end
chx = str2double(chx_string);
if (chx > nch || chx <= 0)
    errordlg('please input correct channel number','parameter lost');
    return
end

% channel Y
chy_string = get(handles.chy,'String');
if chy_string == ' '
    errordlg('please input channel number','parameter lost');
    return
end
chy = str2double(chy_string);
if (chy > nch || chy <= 0)
    errordlg('please input correct channel number','parameter lost');
    return
end
if (chx == chy)
    errordlg('please input different channel number','parameter lost');
    return
end

% core: coherence view
co_view(coh,chx,chy,fs);

% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 delete(handles.coh_view);

% --- Executes on button press in pot_spect.
function pot_spect_Callback(hObject, eventdata, handles)
% hObject    handle to pot_spect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% coherence data
if isequal(handles.filepath,'')
    errordlg('File not found','File Error');
    return
end
data   = load(handles.filepath);
name   = fieldnames(data);
d      = name{1};
coh    = getfield(data,d);
% channel number
c = size(coh,3);
nch = (1+sqrt(1+8*c))/2;

% sampling rate
fs_string = get(handles.fs,'String');
if fs_string == ' '
    errordlg('please input sampling rate','parameter lost');
    return
end
fs = str2double(fs_string);
if (fs <= 0)
    errordlg('please input correct sampling rate','parameter lost');
    return
end

% channel X
chx_string = get(handles.chx,'String');
if chx_string == ' '
    errordlg('please input channel number','parameter lost');
    return
end
chx = str2double(chx_string);
if (chx > nch || chx <= 0)
    errordlg('please input correct channel number','parameter lost');
    return
end

% channel Y
chy_string = get(handles.chy,'String');
if chy_string == ' '
    errordlg('please input channel number','parameter lost');
    return
end
chy = str2double(chy_string);
if (chy > nch || chy <= 0)
    errordlg('please input correct channel number','parameter lost');
    return
end
if (chx == chy)
    errordlg('please input different channel number','parameter lost');
    return
end

% window
win_string=get(handles.win,'String');
if win_string == ' '
    errordlg('please input window','parameter lost');
    return
end
win = str2num(win_string);
% TODO: window error check

% core: coherence view
co_view(coh,chx,chy,fs,win);


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


% --- Executes on button press in choose_coh.
function choose_coh_Callback(hObject, eventdata, handles)
% hObject    handle to choose_coh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile( ...
    {'*.mat','MAT-files (*.mat)'; ...
    '*.*',  'All Files (*.*)'},...
    'Pick a file');

if isequal(filename,0) || isequal(pathname,0)
    filepath = handles.filepath;
    if isempty(filepath)
        filename = 'No file selected...';
    else
        [path,name,ext] = fileparts(filepath);
        filename = [name,ext];
    end%if
else
    handles.filepath = fullfile(pathname,filename);
end

guidata(hObject,handles);
set(handles.fcoh_show,'String',filename);


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


