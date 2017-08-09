function varargout = pairwise_coherence(varargin)
% PAIRWISE_COHERENCE M-file for pairwise_coherence.fig
%      PAIRWISE_COHERENCE, by itself, creates a new PAIRWISE_COHERENCE or raises the existing
%      singleton*.
%
%      H = PAIRWISE_COHERENCE returns the handle to a new PAIRWISE_COHERENCE or the handle to
%      the existing singleton*.
%
%      PAIRWISE_COHERENCE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PAIRWISE_COHERENCE.M with the given input arguments.
%
%      PAIRWISE_COHERENCE('Property','Value',...) creates a new PAIRWISE_COHERENCE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pairwise_coherence_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pairwise_coherence_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.2$ $Date: 12-Sep-2007 15:33:18$
% SHIS UT-Houston, Houston, TX 77030, USA.
% 
% Lei Xu, Hualou Liang

% Edit the above text to modify the response to help pairwise_coherence

% Last Modified by GUIDE v2.5 15-Sep-2007 16:36:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pairwise_coherence_OpeningFcn, ...
                   'gui_OutputFcn',  @pairwise_coherence_OutputFcn, ...
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


% --- Executes just before pairwise_coherence is made visible.
function pairwise_coherence_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pairwise_coherence (see VARARGIN)

% Choose default command line output for pairwise_coherence
handles.output = hObject;

% Update handles structure
handles.filepath = '';
guidata(hObject, handles);

% UIWAIT makes pairwise_coherence wait for user response (see UIRESUME)
% uiwait(handles.pair_coh);


% --- Outputs from this function are returned to the command line.
function varargout = pairwise_coherence_OutputFcn(hObject, eventdata, handles) 
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
else
    handles.filepath=dname;
end

guidata(hObject,handles);
set(handles.dir_show,'String',handles.filepath);

% restore directory
cd(cdir);


% --- Executes on button press in coh_bi.
function coh_bi_Callback(hObject, eventdata, handles)
% hObject    handle to coh_bi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(handles.filepath,'')
    errordlg('Directory not found','Directory Error');
    return
end

fs_string=get(handles.fs,'String');
n_string=get(handles.nbin,'String');   
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

% core: coherence bivariate model
coh = bi_coherence(fpath,nbin,fs); %#ok<NASGU>

% save
suggestname = 'coh_bi';
[filename,pathname] = uiputfile( ...
    {'*.mat','MAT-files (*.mat)';...
    '*.*',  'All Files (*.*)'},...
    'Save as',...
    suggestname);
if isequal(filename,0) || isequal(pathname,0)
else
    filepath = fullfile(pathname,filename);
    save(filepath,'coh');
end%if


% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 delete(handles.pair_coh);


function dir_show_Callback(hObject, eventdata, handles)
% hObject    handle to dir_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dir_show as text
%        str2double(get(hObject,'String')) returns contents of dir_show as a double


% --- Executes during object creation, after setting all properties.
function dir_show_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dir_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




