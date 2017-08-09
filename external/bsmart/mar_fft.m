function varargout = mar_fft(varargin)
% MAR_FFT M-file for mar_fft.fig
%      MAR_FFT, by itself, creates a new MAR_FFT or raises the existing
%      singleton*.
%
%      H = MAR_FFT returns the handle to a new MAR_FFT or the handle to
%      the existing singleton*.
%
%      MAR_FFT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAR_FFT.M with the given input arguments.
%
%      MAR_FFT('Property','Value',...) creates a new MAR_FFT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mar_fft_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mar_fft_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mar_fft

% Last Modified by GUIDE v2.5 14-Sep-2007 15:59:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mar_fft_OpeningFcn, ...
                   'gui_OutputFcn',  @mar_fft_OutputFcn, ...
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


% --- Executes just before mar_fft is made visible.
function mar_fft_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mar_fft (see VARARGIN)

% Choose default command line output for mar_fft
handles.output = hObject;

% Update handles structure
handles.nchannels=0;
handles.ntrails=0;
handles.npoints=0;
handles.channel=0;
handles.trail=0;
handles.n=0;
handles.fs=0;
handles.filepath='';
guidata(hObject, handles);

% UIWAIT makes mar_fft wait for user response (see UIRESUME)
% uiwait(handles.FFT);


% --- Outputs from this function are returned to the command line.
function varargout = mar_fft_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
{ '*.mat','MAT-files (*.mat)'; ...
'*.*',  'All Files (*.*)'; ...
'*.m;*.fig;*.mat;*.mdl','MATLAB Files (*.m,*.fig,*.mat,*.mdl)';...
   '*.m',  'M-files (*.m)'; ...
   '*.fig','Figures (*.fig)'; ...
   '*.mdl','Models (*.mdl)'; ...
   '*.dm6','steven file'},...
   'Pick a file');
if isequal(filename,0) | isequal(pathname,0)
else
handles.filepath=fullfile(pathname,filename);
end
guidata(hObject,handles);
set(handles.edit1,'String',handles.filepath);


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fs_string=get(handles.edit11,'String');
n_string=get(handles.edit9,'String');
channelchoose_string=get(handles.edit8,'String');
trailchoose_string=get(handles.edit10,'String');
trailchoose_string2=get(handles.edit14,'String');
pointchoose_string=get(handles.edit15,'String');
pointchoose_string2=get(handles.edit16,'String');
if isequal(handles.filepath,'')
errordlg('File not found','File Error');
return;
end
data=load(handles.filepath);
name=fieldnames(data);
d=name{1,1};
dat=getfield(data,d);
si=size(dat);
if fs_string==' '
       errordlg('please input sampling frequency','parameter lost');return; end
   if channelchoose_string==' '
       errordlg('please input the channel','parameter lost');return; end
   if trailchoose_string==' '
       errordlg('please input the start trial','parameter lost');return; end
   if trailchoose_string2==' '
       errordlg('please input the end trial','parameter lost');return; end
   if pointchoose_string==' '
       errordlg('please input the start point','parameter lost');return; end
   if pointchoose_string2==' '
       errordlg('please input the end point','parameter lost');return; end
   if n_string==' '
       errordlg('please input the length of fft','parameter lost');return; end
   handles.fs=str2num(get(handles.edit11,'String'));
   handles.npoints=si(1);
   handles.ntrails=si(3);
   handles.nchannels=si(2);
   handles.channel=str2num(get(handles.edit8,'String'));
   handles.trail=str2num(get(handles.edit10,'String'));
   handles.trail2=str2num(get(handles.edit14,'String'));
   handles.point=str2num(get(handles.edit15,'String'));
   handles.point2=str2num(get(handles.edit16,'String'));
   handles.n=str2num(get(handles.edit9,'String'));
   guidata(hObject,handles);
   
   if handles.channel>handles.nchannels || handles.channel<=0
       errordlg('please input correct channel','parameter lost');return; end
   if handles.point>handles.npoints || handles.point<=0 || handles.point2<handles.point || handles.point2<0 || handles.point2>handles.npoints
       errordlg('please input correct point','parameter lost');return; end
   if handles.trail>handles.ntrails || handles.trail<=0 || handles.trail2<handles.trail || handles.trail2<0 || handles.trail2>handles.ntrails
       errordlg('please input correct trial','parameter lost');return; end
   fs=handles.fs;
   channel=handles.nchannels;
   choosechannel=handles.channel;
   trail=handles.ntrails;
   choosetrail=handles.trail;
   choosetrail2=handles.trail2;
   choosepoint=handles.point;
   choosepoint2=handles.point2;
   points=handles.npoints;
   n=handles.n;
   %a=dat(:,choosechannel,choosetrail:choosetrail2);
   yy=[];
    figure('Name','FFT','NumberTitle','off')
   for i=choosetrail:choosetrail2
   y=fft(dat(choosepoint:choosepoint2,choosechannel,i),n);
   pp=20*log10(y.*conj(y)/n);
   size(pp);
   yy=[yy pp];
   f=fs*(0:n/2)/n;
    plot(f,pp(1:n/2+1));
    hold on
   end
   yy=mean(yy,2);
   plot(fs*(0:n/2)/n,yy(1:n/2+1),'r');
   h=gca;
   xlabel(h,'frequency (Hz)');
   ylabel(h,'Power Spectrum (db)');
   title('Red curve is the average');
   figure('Name','Average_FFT','NumberTitle','off')
   plot(fs*(0:n/2)/n,yy(1:n/2+1),'r');
   h=gca;
   xlabel(h,'frequency (Hz)');
   ylabel(h,'Power Spectrum (db)');
  
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
                delete(handles.FFT);
function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


