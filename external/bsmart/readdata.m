function varargout = readdata(varargin)
% READDATA M-file for readdata.fig
%      READDATA, by itself, creates a new READDATA or raises the existing
%      singleton*.
%
%      H = READDATA returns the handle to a new READDATA or the handle to
%      the existing singleton*.
%
%      READDATA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in READDATA.M with the given input arguments.
%
%      READDATA('Property','Value',...) creates a new READDATA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before readdata_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to readdata_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help readdata

% Last Modified by GUIDE v2.5 14-Sep-2007 14:32:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @readdata_OpeningFcn, ...
                   'gui_OutputFcn',  @readdata_OutputFcn, ...
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


% --- Executes just before readdata is made visible.
function readdata_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to readdata (see VARARGIN)

% Choose default command line output for readdata
handles.output = hObject;

% Update handles structure
handles.nchannels=0;
handles.ntrails=0;
handles.npoints=0;
handles.filepath='';
guidata(hObject, handles);

% UIWAIT makes readdata wait for user response (see UIRESUME)
% uiwait(handles.import);


% --- Outputs from this function are returned to the command line.
function varargout = readdata_OutputFcn(hObject, eventdata, handles)  %#ok<INUSL>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in choose_file.
function choose_file_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to choose_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% select a file
[filename,pathname] = uigetfile( ...
    {'*.bin','Binary Files (*.bin)';...
    '*.dm6','Steven Files (*.dm6)';...
    '*.*',  'All Files (*.*)'},...
    'Pick a file');

if isequal(filename,0) || isequal(pathname,0)
    filename = 'No file selected';
else
    handles.filepath = fullfile(pathname,filename);
end

% show file name
set(handles.fname_show,'String',filename);

% save data into handles
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function fname_show_CreateFcn(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
% hObject    handle to fname_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.


function input_chan_Callback(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
% hObject    handle to input_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_chan as text
%        str2double(get(hObject,'String')) returns contents of input_chan as a double


% --- Executes during object creation, after setting all properties.
function input_chan_CreateFcn(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
% hObject    handle to input_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function input_trl_Callback(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
% hObject    handle to input_trl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_trl as text
%        str2double(get(hObject,'String')) returns contents of input_trl as a double


% --- Executes during object creation, after setting all properties.
function input_trl_CreateFcn(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
% hObject    handle to input_trl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function input_pts_Callback(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
% hObject    handle to input_pts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_pts as text
%        str2double(get(hObject,'String')) returns contents of input_pts as a double


% --- Executes during object creation, after setting all properties.
function input_pts_CreateFcn(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
% hObject    handle to input_pts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in import_data.
function import_data_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to import_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% check error
if isequal(handles.filepath,'')
    errordlg('File not found','File Error');
    return
end

point_string   = get(handles.input_pts,'String');
channel_string = get(handles.input_chan,'String');
trial_string   = get(handles.input_trl,'String');
    
% error check
if point_string == ' '
    errordlg('Please input number of points','Parameter lost');
    return
end
if channel_string == ' '
    errordlg('Please input number of channels','Parameter lost');
    return
end
if  trial_string == ' '
    errordlg('Please input number of trials','Parameter lost');
    return
end

fpath  = handles.filepath;
npts   = str2double(point_string);
nchn   = str2double(channel_string);
ntrl   = str2double(trial_string);

% core: data format conversion
dat = readdat(fpath,npts,nchn,ntrl); %#ok<NASGU>

% save file
oldfullname = handles.filepath;
[pathstr, suggestname] = fileparts(oldfullname);
[filename,pathname] = uiputfile( ...
    {'*.mat','MAT-files (*.mat)';...
    '*.*',  'All Files (*.*)'},...
    'Save as',...
    suggestname);

if isequal(filename,0) || isequal(pathname,0)
else
    filepath = fullfile(pathname,filename);
    save(filepath,'dat');
end%if

% save data in handles structure
handles.npoints   = npts;
handles.nchaanels = nchn;
handles.ntrials   = ntrl;
guidata(hObject,handles);


% --- Executes on button press in close_import.
function close_import_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to close_import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.import);

% [EOF]
