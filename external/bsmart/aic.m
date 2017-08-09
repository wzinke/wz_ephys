function varargout = aic(varargin)
% AIC M-file for aic.fig
%      AIC, by itself, creates a new AIC or raises the existing
%      singleton*.
%
%      H = AIC returns the handle to a new AIC or the handle to
%      the existing singleton*.
%
%      AIC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AIC.M with the given input arguments.
%
%      AIC('Property','Value',...) creates a new AIC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before aic_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to aic_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

%   Copyright 2007 BSMART group
%   $Revision: 0.1 $  $Date: 28-Aug-2007 10:36:41 $
%   by Richard J. Cui
%   SHIS UT-Houston, Houston, TX 77030, USA.

% Edit the above text to modify the response to help aic

% Last Modified by GUIDE v2.5 14-Sep-2007 16:10:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @aic_OpeningFcn, ...
                   'gui_OutputFcn',  @aic_OutputFcn, ...
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


% --- Executes just before aic is made visible.
function aic_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to aic (see VARARGIN)

% Choose default command line output for aic
handles.output = hObject;

% Update handles structure
handles.nchannels=0;
handles.ntrails=0;
handles.npoints=0;
handles.window=0;
handles.order=0;
handles.filepath='';
guidata(hObject, handles);

% UIWAIT makes aic wait for user response (see UIRESUME)
% uiwait(handles.AIC);


% --- Outputs from this function are returned to the command line.
function varargout = aic_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function choose_mat_file_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile( ...
    { '*.mat','MAT-files (*.mat)'; ...
    '*.*',  'All Files (*.*)'},...
    'Pick a file');
if isequal(filename,0) || isequal(pathname,0)
else
    handles.filepath=fullfile(pathname,filename);
end
guidata(hObject,handles);
set(handles.file_show,'String',handles.filepath);


function file_show_Callback(hObject, eventdata, handles)
% hObject    handle to file_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of file_show as text
%        str2double(get(hObject,'String')) returns contents of file_show as a double


% --- Executes during object creation, after setting all properties.
function file_show_CreateFcn(hObject, eventdata, handles)
% hObject    handle to file_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in aic_test.
function aic_test_Callback(hObject, eventdata, handles)
% hObject    handle to aic_test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

window_string = get(handles.winlen,'String');
order_string  = get(handles.maxorder,'String');
% error check
if isequal(handles.filepath,'')
    errordlg('File not found','File Error');
    return
end

datafile = handles.filepath;
data     = load(datafile);
name     = fieldnames(data);
d        = name{1};

% core: aic_test
dat      = getfield(data,d); %#ok<GFLD>
winlen   = str2double(window_string);
maxorder = str2double(order_string);
points   = size(dat,1);
aic = aic_test(dat,winlen,maxorder);

% save AIC
suggestname = 'AIC';
[filename,pathname] = uiputfile( ...
    {'*.mat','MAT-files (*.mat)';...
    '*.*',  'All Files (*.*)'},...
    'Save as',...
    suggestname);

if isequal(filename,0) || isequal(pathname,0)
else
    filepath = fullfile(pathname,filename);
    save(filepath,'aic');
end%if

AIC_order_view(points,winlen,aic);

   
% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.AIC);


function winlen_Callback(hObject, eventdata, handles)
% hObject    handle to winlen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of winlen as text
%        str2double(get(hObject,'String')) returns contents of winlen as a double


% --- Executes during object creation, after setting all properties.
function winlen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to winlen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function maxorder_Callback(hObject, eventdata, handles)
% hObject    handle to maxorder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxorder as text
%        str2double(get(hObject,'String')) returns contents of maxorder as a double


% --- Executes during object creation, after setting all properties.
function maxorder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxorder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

