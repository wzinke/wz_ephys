function varargout = Lyapunov(varargin)
% LYAPUNOV M-file for Lyapunov.fig
%      LYAPUNOV, by itself, creates a new LYAPUNOV or raises the existing
%      singleton*.
%
%      H = LYAPUNOV returns the handle to a new LYAPUNOV or the handle to
%      the existing singleton*.
%
%      LYAPUNOV('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LYAPUNOV.M with the given input arguments.
%
%      LYAPUNOV('Property','Value',...) creates a new LYAPUNOV or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Lyapunov_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Lyapunov_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help Lyapunov

% Last Modified by GUIDE v2.5 15-Sep-2007 12:36:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Lyapunov_OpeningFcn, ...
                   'gui_OutputFcn',  @Lyapunov_OutputFcn, ...
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


% --- Executes just before Lyapunov is made visible.
function Lyapunov_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Lyapunov (see VARARGIN)

% Choose default command line output for Lyapunov
handles.output = hObject;

% Update handles structure
handles.arcfilepath = '';
handles.arnfilepath = '';
guidata(hObject, handles);

% UIWAIT makes Lyapunov wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Lyapunov_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in choose_farc_show.
function choose_arc_Callback(hObject, eventdata, handles)
% hObject    handle to choose_farc_show (see GCBO)
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
    handles.arcfilepath = fullfile(pathname,filename);
end
guidata(hObject,handles);
set(handles.farc_show,'String',handles.arcfilepath);

% --- Executes on button press in choose_nvm.
function choose_nvm_Callback(hObject, eventdata, handles)
% hObject    handle to choose_nvm (see GCBO)
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



% --- Executes on button press in sta_test.
function sta_test_Callback(hObject, eventdata, handles)
% hObject    handle to sta_test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(handles.arcfilepath,'')
    errordlg('Coffecient File not found','File Error');
    return
end
if isequal(handles.arnfilepath,'')
    errordlg('Noise File not found','File Error');
    return
end

t_string = get(handles.Tgen,'String');
if t_string == ' '
    errordlg('Please input T (points to be generated)','parameter lost');
    return
end

% AR coefficients
a    = load(handles.arcfilepath);
name = fieldnames(a);
d    = name{1};
arc  = getfield(a,d); %#ok<GFLD>
% make it compatable with Matlab R14
if size(arc,2) == 1
    arc = arc.';    % if arc is a column, change it to row
end%if

% Noise cov mat
a    = load(handles.arnfilepath);
name = fieldnames(a);
d    = name{1};
ncm  = getfield(a,d); %#ok<GFLD>
% make it compatable with Matlab R14
if size(ncm,2) == 1
    ncm = ncm.';    % if ncm is a column, change it to row
end%if

% points to be generated
T = str2double(t_string);
le = lyap_batch(arc,ncm,T);

% save results
suggestname = 'Lyap';
[filename,pathname] = uiputfile( ...
    {'*.mat','MAT-files (*.mat)';...
    '*.*',  'All Files (*.*)'},...
    'Save as',...
    suggestname);
if isequal(filename,0) || isequal(pathname,0)
else
    filepath = fullfile(pathname,filename);
    save(filepath,'le');
end%if

% plot results
figure('Name','Stability Test','NumberTitle','off')
plot(le,'-o');
hold on;
ax = axis;
plot([ax(1) ax(2)],[0 0],':');
h = gca;
xlabel(h,'Time')
ylabel(h,'Stability Index')
set(h,'YLim',[-0.2 0.2]);


% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.Lyapunov);



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



function Tgen_Callback(hObject, eventdata, handles)
% hObject    handle to Tgen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tgen as text
%        str2double(get(hObject,'String')) returns contents of Tgen as a double


% --- Executes during object creation, after setting all properties.
function Tgen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tgen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


