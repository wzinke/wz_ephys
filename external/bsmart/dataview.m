function varargout = dataview(varargin)
% DATAVIEW M-file for dataview.fig
%      DATAVIEW, by itself, creates a new DATAVIEW or raises the existing
%      singleton*.
%
%      H = DATAVIEW returns the handle to a new DATAVIEW or the handle to
%      the existing singleton*.
%
%      DATAVIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DATAVIEW.M with the given input arguments.
%
%      DATAVIEW('Property','Value',...) creates a new DATAVIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dataview_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dataview_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu_file.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

%   Copyright 2007 BSMART group.
%   $Revision: 0.2 $  $Date: 15-Sep-2007 22:43:27 $
%   SHIS UT-Houston, Houston, TX 77030, USA.

% Edit the above text to modify the response to help dataview

% Last Modified by GUIDE v2.5 15-Sep-2007 22:39:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dataview_OpeningFcn, ...
                   'gui_OutputFcn',  @dataview_OutputFcn, ...
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


% --- Executes just before dataview is made visible.
function dataview_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dataview (see VARARGIN)

% Choose default command line output for dataview
handles.output = hObject;

% Update handles structure
handles.nchannels=0;
handles.ntrails=0;
handles.npoints=0;
handles.filepath='';
guidata(hObject, handles);

% UIWAIT makes dataview wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dataview_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function fdata_show_Callback(hObject, eventdata, handles)
% hObject    handle to fdata_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fdata_show as text
%        str2double(get(hObject,'String')) returns contents of fdata_show as a double


% --- Executes during object creation, after setting all properties.
function fdata_show_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fdata_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in chartview.
function chartview_Callback(hObject, eventdata, handles)
% hObject    handle to chartview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% check file validation
if isequal(handles.filepath,'')
    errordlg('File not found','File Error');
    return;
end
% get data
data2 = load(handles.filepath);
name = fieldnames(data2);
d = name{1,1};
dat = getfield(data2,d); %#ok<GFLD>
si = size(dat);
handles.nchannels = si(2);
handles.ntrails = si(3);
handles.npoints = si(1);
guidata(hObject,handles);
% global channel trail point data;
data = dat;
% channel = handles.nchannels;
% trail = handles.ntrails;
% point =handles.npoints;
% chart_view;
% change data format to channel x timepoints x trials
sigdata = permute(data,[2 1 3]);
sigplot(sigdata);

% --- Executes on button press in grid_view.
function grid_view_Callback(hObject, eventdata, handles)
% hObject    handle to grid_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(handles.filepath,'')
    errordlg('File not found','File Error');
    return
end

data2=load(handles.filepath);
name=fieldnames(data2);
d=name{1,1};
dat=getfield(data2,d);

% 
%     si=size(dat);
%     handles.nchannels=si(2);
%     handles.ntrails=si(3);
%     handles.npoints=si(1);
% guidata(hObject,handles);
% global channel trail point data;
% channel  = handles.nchannels;
% trail    = handles.ntrails;
% point    = handles.npoints;

grid_view(dat);


% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);



% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function file_open_Callback(hObject, eventdata, handles)
% hObject    handle to file_open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile( ...
    {'*.mat','MAT-files (*.mat)'; ...
    '*.*',  'All Files (*.*)'},...
    'Pick a file');
if isequal(filename,0) || isequal(pathname,0)
    filepath = handles.filepath;
    if isempty(filepath)
        filename = 'No data file selected...';
    else
        [pathstr,name,ext,versn] = fileparts(filepath); %#ok<ASGLU,NASGU>
        filename = [name,ext];
    end%if
else
    handles.filepath = fullfile(pathname,filename);
end
guidata(hObject,handles);

set(handles.fdata_show,'String',filename);



% --------------------------------------------------------------------
function menu_help_Callback(hObject, eventdata, handles)
% hObject    handle to menu_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in button_open.
function button_open_Callback(hObject, eventdata, handles)
% hObject    handle to button_open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

file_open_Callback(hObject, eventdata, handles);
