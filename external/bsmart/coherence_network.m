function varargout = coherence_network(varargin)
% COHERENCE_NETWORK M-file for coherence_network.fig
%      COHERENCE_NETWORK, by itself, creates a new COHERENCE_NETWORK or raises the existing
%      singleton*.
%
%      H = COHERENCE_NETWORK returns the handle to a new COHERENCE_NETWORK or the handle to
%      the existing singleton*.
%
%      COHERENCE_NETWORK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COHERENCE_NETWORK.M with the given input arguments.
%
%      COHERENCE_NETWORK('Property','Value',...) creates a new COHERENCE_NETWORK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before coherence_network_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to coherence_network_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help coherence_network

% Last Modified by GUIDE v2.5 15-Sep-2007 18:44:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @coherence_network_OpeningFcn, ...
                   'gui_OutputFcn',  @coherence_network_OutputFcn, ...
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


% --- Executes just before coherence_network is made visible.
function coherence_network_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to coherence_network (see VARARGIN)

% Choose default command line output for coherence_network
handles.output = hObject;

% Update handles structure
handles.filepath='';
handles.filepath2='';
guidata(hObject, handles);

% UIWAIT makes coherence_network wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = coherence_network_OutputFcn(hObject, eventdata, handles) 
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
if isequal(filename,0) | isequal(pathname,0)
else
handles.filepath=fullfile(pathname,filename);end
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


function startf_Callback(hObject, eventdata, handles)
% hObject    handle to startf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startf as text
%        str2double(get(hObject,'String')) returns contents of startf as a double


% --- Executes during object creation, after setting all properties.
function startf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function endf_Callback(hObject, eventdata, handles)
% hObject    handle to endf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endf as text
%        str2double(get(hObject,'String')) returns contents of endf as a double


% --- Executes during object creation, after setting all properties.
function endf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to endf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function window_Callback(hObject, eventdata, handles)
% hObject    handle to window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of window as text
%        str2double(get(hObject,'String')) returns contents of window as a double


% --- Executes during object creation, after setting all properties.
function window_CreateFcn(hObject, eventdata, handles)
% hObject    handle to window (see GCBO)
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

delete(handles.figure1);

% --- Executes on button press in show_net.
function show_net_Callback(hObject, eventdata, handles)
% hObject    handle to show_net (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% coherence file
if isequal(handles.filepath,'')
    errordlg('File not found','File Error');
    return
end
a    = load(handles.filepath);
name = fieldnames(a);
d    = name{1};
coh  = getfield(a,d);

% location file
if isequal(handles.filepath2,'')
    errordlg('File not found','File Error');
    return
end
a    = load(handles.filepath2);
name = fieldnames(a);
d    = name{1};
loc  = getfield(a,d);

fre1_string=get(handles.startf,'String');
if fre1_string==' '
    errordlg('please input correct frequency','parameter lost');
    return
end
fstart = str2double(fre1_string);
if (fstart <= 0)
    errordlg('please input correct frequency','parameter lost');
    return
end

fre2_string=get(handles.endf,'String');
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

time_string = get(handles.window,'String');
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
s1  = size(coh);
c   = s1(3);
nch = (1+sqrt(1+8*c))/2;    % number of channels
chint_string = get(handles.chan_int,'String');
if isempty(chint_string)
    chan = 1:nch;
else
    chan = str2num(chint_string); %#ok<ST2NM>
end

% core: shaw coh network
conetwork(coh,loc,thre,win,fstart,fend,chan)


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


% --- Executes on button press in loc_mat_file.
function loc_mat_file_Callback(hObject, eventdata, handles)
% hObject    handle to loc_mat_file (see GCBO)
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
handles.filepath2=fullfile(pathname,filename);end
guidata(hObject,handles);
set(handles.floc_show,'String',handles.filepath2);


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



function chan_int_Callback(hObject, eventdata, handles)
% hObject    handle to chan_int (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of chan_int as text
%        str2double(get(hObject,'String')) returns contents of chan_int as a double


% --- Executes during object creation, after setting all properties.
function chan_int_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chan_int (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


