function varargout = power_pairwise(varargin)
% POWER_PAIRWISE M-file for power_pairwise.fig
%      POWER_PAIRWISE, by itself, creates a new POWER_PAIRWISE or raises the existing
%      singleton*.
%
%      H = POWER_PAIRWISE returns the handle to a new POWER_PAIRWISE or the handle to
%      the existing singleton*.
%
%      POWER_PAIRWISE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POWER_PAIRWISE.M with the given input arguments.
%
%      POWER_PAIRWISE('Property','Value',...) creates a new POWER_PAIRWISE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before power_pairwise_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to power_pairwise_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.3$ $Date: 13-Sep-2007 17:04:26$
% SHIS UT-Houston, Houston, TX 77030, USA.
% 
% Lei Xu, Hualou Liang


% Edit the above text to modify the response to help power_pairwise

% Last Modified by GUIDE v2.5 15-Sep-2007 12:56:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @power_pairwise_OpeningFcn, ...
                   'gui_OutputFcn',  @power_pairwise_OutputFcn, ...
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


% --- Executes just before power_pairwise is made visible.
function power_pairwise_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to power_pairwise (see VARARGIN)

% Choose default command line output for power_pairwise
handles.output = hObject;

% Update handles structure
handles.filepath = '';
guidata(hObject, handles);

% UIWAIT makes power_pairwise wait for user response (see UIRESUME)
% uiwait(handles.power_pairwise);


% --- Outputs from this function are returned to the command line.
function varargout = power_pairwise_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



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



function bin_Callback(hObject, eventdata, handles)
% hObject    handle to bin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bin as text
%        str2double(get(hObject,'String')) returns contents of bin as a double


% --- Executes during object creation, after setting all properties.
function bin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in power.
function power_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(handles.filepath,'')
    errordlg('Directory not found','Directory Error');
    return;
end

fs_string = get(handles.fs,'String');
n_string  = get(handles.bin,'String');   
if fs_string == ' '
         errordlg('please input sample frequency','parameter lost');
         return
end
if n_string == ' '
         errordlg('please input the length of returned spectrum','parameter lost');
         return
end

fpath = handles.filepath;
nbin  = str2double(n_string);
fs    = str2double(fs_string);

% core: autopower from bivariate model
power = bi_power(fpath,nbin,fs); %#ok<NASGU>

% save power
suggestname = 'power_bi';
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

 delete(handles.power_pairwise);

function path_show_Callback(hObject, eventdata, handles)
% hObject    handle to path_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of path_show as text
%        str2double(get(hObject,'String')) returns contents of path_show as a double


% --- Executes during object creation, after setting all properties.
function path_show_CreateFcn(hObject, eventdata, handles)
% hObject    handle to path_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --- Executes on button press in choose_dir.
function choose_dir_Callback(hObject, eventdata, handles)
% hObject    handle to choose_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% change directory
cdir = pwd;     % current directory
path = mfilename('fullpath');
fdir = fileparts(path);

dname = uigetdir(fdir,'Select a directory');
if isequal(dname,0)
    dname = handles.filepath;
    if isempty(dname)
        dname = 'No directory selected...';
    end%if
else
    handles.filepath = dname;
end
guidata(hObject,handles);
set(handles.path_show,'String',dname);

% restore directory
cd(cdir);
