function varargout = moving_window_pairwise(varargin)
% MOVING_WINDOW_PAIRWISE M-file for moving_window_pairwise.fig
%      MOVING_WINDOW_PAIRWISE, by itself, creates a new MOVING_WINDOW_PAIRWISE or raises the existing
%      singleton*.
%
%      H = MOVING_WINDOW_PAIRWISE returns the handle to a new MOVING_WINDOW_PAIRWISE or the handle to
%      the existing singleton*.
%
%      MOVING_WINDOW_PAIRWISE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MOVING_WINDOW_PAIRWISE.M with the given input arguments.
%
%      MOVING_WINDOW_PAIRWISE('Property','Value',...) creates a new MOVING_WINDOW_PAIRWISE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before moving_window_pairwise_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to moving_window_pairwise_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright (c) 2006-2007 BSMART Group
% by Richard Cui
% $Revision: 0.2$ $Date: 14-Sep-2007 20:27:40$
% SHIS UT-Houston, Houston, TX 77030, USA.
% 
% Lei Xu, Hualou Liang


% Edit the above text to modify the response to help moving_window_pairwise

% Last Modified by GUIDE v2.5 14-Sep-2007 20:20:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @moving_window_pairwise_OpeningFcn, ...
                   'gui_OutputFcn',  @moving_window_pairwise_OutputFcn, ...
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


% --- Executes just before moving_window_pairwise is made visible.
function moving_window_pairwise_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to moving_window_pairwise (see VARARGIN)

% Choose default command line output for moving_window_pairwise
handles.output = hObject;

% Update handles structure
handles.nchannels=0;
handles.ntrails=0;
handles.npoints=0;
handles.window=0;
handles.order=0;
handles.filepath='';
% handles.fs=0;
% handles.n=0;
% handles.freq=[];
guidata(hObject, handles);

% UIWAIT makes moving_window_pairwise wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = moving_window_pairwise_OutputFcn(hObject, eventdata, handles) 
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
'*.m;*.fig;*.mat;*.mdl','MATLAB Files (*.m,*.fig,*.mat,*.mdl)';
   '*.m',  'M-files (*.m)'; ...
   '*.fig','Figures (*.fig)'; ...
   '*.mdl','Models (*.mdl)'; ...
   '*.*',  'All Files (*.*)'; ...
   '*.dm6','steven file'},...
   'Pick a file');
if isequal(filename,0) | isequal(pathname,0)
   %disp('User selected Cancel')
else
handles.filepath=fullfile(pathname,filename);
end
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
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function morder_Callback(hObject, eventdata, handles)
% hObject    handle to morder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of morder as text
%        str2double(get(hObject,'String')) returns contents of morder as a double


% --- Executes during object creation, after setting all properties.
function morder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to morder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% start_string=get(handles.spts,'String');
% end_string=get(handles.epts,'String');
% order_string=get(handles.morder,'String');
% fs_string=get(handles.edit3,'String');
% window_string=get(handles.winlen,'String');
% n_string=get(handles.edit4,'String');
% if isequal(handles.filepath,'')
% errordlg('File not found','File Error');
% return;
% end
% data=load(handles.filepath);
% name=fieldnames(data);
% d=name{1,1};
% dat=getfield(data,d);
% si=size(dat);
%    if start_string==' '
%         errordlg('please input the starting position','parameter lost');return; end
%    if end_string==' '
%         errordlg('please input the ending position','parameter lost');return; end
%    if window_string==' '
%         errordlg('please input width of window','parameter lost');return; end
%    if order_string==' '
%         errordlg('please input the model order','parameter lost');return; end
%    if fs_string==' '
%         errordlg('please input sample frequency','parameter lost');return; end
%    if n_string==' '
%         errordlg('please input the length of returned spectrum','parameter lost');return; end
%    handles.start=str2num(get(handles.spts,'String'));
%    handles.end=str2num(get(handles.epts,'String'));
%    handles.fs=str2num(get(handles.edit3,'String'));
%    handles.n=str2num(get(handles.edit4,'String'));
%    handles.window=str2num(get(handles.winlen,'String'));
%    handles.npoints=si(1);
%    handles.ntrails=si(3);
%    handles.nchannels=si(2);
%    handles.start=str2num(get(handles.spts,'String'));
%    handles.end=str2num(get(handles.epts,'String'));
%    handles.order=str2num(get(handles.morder,'String'));
%    guidata(hObject,handles);
%    n=handles.n;
%    fs=handles.fs;
%    if handles.start>handles.npoints
%        errordlg('please input the correct starting position','parameter lost');return; end
%    if handles.end<handles.start || handles.end>handles.npoints
%        errordlg('please input the correct ending position','parameter lost');return; end
%    if handles.window>handles.end-handles.start+1
%        errordlg('please input the correct window length','parameter lost');return; end
%    if handles.order>20
%        errordlg('please input the model order less than 20','parameter lost');return; end
%    channel2=handles.nchannels;
%    
%    channel=2;
%    save channel channel -ascii;
%    
%    trail=handles.ntrails;
%    save trail trail -ascii;
%    
%    points=handles.end-handles.start+1;
%    save points points -ascii;
%    
%    window=handles.window;
%    save window window -ascii;
%    
%    order=handles.order;
%    save order order -ascii;
%    dat=dat(handles.start:handles.end,:,:);
%    %pairco=zeros(points-window+1,n,channel2*(channel2-1)/2);
%    %spect=zeros(points-window+1,n,channel2);
%    pairco=[];
%    spect=[];
%    k=0;l=1;
%    for i=1:(channel2-1)
%        for j=(i+1):channel2
%            dat1=dat(:,i,:);
%            dat2=dat(:,j,:);
%            dat3=cat(2,dat1,dat2);
%            writedat('dataset.bin',dat3);
%            eval(['unix ' '(''' 'opssmov ' 'dataset.bin ' ' A ' 'Ve' ''')']);
%            aredat=load('A');
%            arndat=load('Ve');
%            [csd4d, freq] = MAR_csd4dmem(aredat,arndat,n,fs);
%            [paircoh, partcoh, mulcoh, autospect] = MAR_coh3D(csd4d, fs);
%            k=k+1;
%            pairco(:,:,k)=paircoh(:,:,1);
% 
%            if i==1
%                l=l+1;
%            spect(:,:,i)=autospect(:,:,1);
%            spect(:,:,l)=autospect(:,:,2);
%     
%            end
%            !DEL A
%            !DEL Ve
%            !DEL dataset.bin
% 
%        end
% 
%    end
%    uisave pairco pairco;
%    uisave spect power;
%     !DEL channel
%     !DEL trail
%     !DEL points
%     !DEL order
%     !DEL window
% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.figure1);

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
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
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
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global Fxy Fyx channel;
% start_string=get(handles.spts,'String');
% end_string=get(handles.epts,'String');
% order_string=get(handles.morder,'String');
% fs_string=get(handles.edit3,'String');
% window_string=get(handles.winlen,'String');
% freq_string=get(handles.edit5,'String');
% if isequal(handles.filepath,'')
% errordlg('File not found','File Error');
% return;
% end
% data=load(handles.filepath);
% name=fieldnames(data);
% d=name{1,1};
% dat=getfield(data,d);
% si=size(dat);
%    if start_string==' '
%         errordlg('please input the starting position','parameter lost');return; end
%    if end_string==' '
%         errordlg('please input the ending position','parameter lost');return; end
%    if window_string==' '
%         errordlg('please input width of window','parameter lost');return; end
%    if order_string==' '
%         errordlg('please input the model order','parameter lost');return; end
%    if fs_string==' '
%         errordlg('please input sample frequency','parameter lost');return; end
%    if freq_string==' '
%         errordlg('please input a vector of frequencies of interest, usually freq=0:fs/2','parameter lost');return; end
%    handles.start=str2num(get(handles.spts,'String'));
%    handles.end=str2num(get(handles.epts,'String'));
%    handles.fs=str2num(get(handles.edit3,'String'));
%    handles.freq=str2num(get(handles.edit5,'String'));
%    handles.window=str2num(get(handles.winlen,'String'));
%    handles.npoints=si(1);
%    handles.ntrails=si(3);
%    handles.nchannels=si(2);
%    handles.order=str2num(get(handles.morder,'String'));
%    guidata(hObject,handles);
%    if handles.start>handles.npoints
%        errordlg('please input the correct starting position','parameter lost');return; end
%    if handles.end<handles.start || handles.end>handles.npoints
%        errordlg('please input the correct ending position','parameter lost');return; end
%    if handles.window>handles.end-handles.start+1
%        errordlg('please input the correct window length','parameter lost');return; end
%    if handles.order>20
%        errordlg('please input the model order less than 20','parameter lost');return; end
%    channel=handles.nchannels;
%    trail=handles.ntrails;
%    points=handles.end-handles.start+1;
%    order=handles.order;
%    fs=handles.fs;
%    freq=handles.freq;
%    win=handles.window;
%    
%    if win>points
%        errordlg('please input the window length no longer than time series','parameter lost');return; end
%    a=zeros(win,channel,trail);
%    b=zeros(channel,trail*win);
%    for i=1:points-win+1
%        a=dat(i+handles.start-1:i+handles.start+win-2,:,:);
%        for k=1:trail
%            for j=1:channel
%                b(j,(k-1)*win+1:k*win)=a(:,j,k);
%            end
%        end
%        [pp,cohe,Fxy(:,:,i),Fyx(:,:,i)]=pwcausal(b,trail,win,order,fs,freq);
%    end
%    uisave Fxy Fxy;
%    uisave Fyx Fyx;
%    granger_causality_movingwindow;
function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



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





function spts_Callback(hObject, eventdata, handles)
% hObject    handle to spts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spts as text
%        str2double(get(hObject,'String')) returns contents of spts as a double


% --- Executes during object creation, after setting all properties.
function spts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function epts_Callback(hObject, eventdata, handles)
% hObject    handle to epts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of epts as text
%        str2double(get(hObject,'String')) returns contents of epts as a double


% --- Executes during object creation, after setting all properties.
function epts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --- Executes on button press in create_model.
function create_model_Callback(hObject, eventdata, handles)
% hObject    handle to create_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(handles.filepath,'')
    errordlg('File not found','File Error');
    return
end

start_string  = get(handles.spts,'String');
end_string    = get(handles.epts,'String');
window_string = get(handles.winlen,'String');
order_string  = get(handles.morder,'String');
if start_string==' '
    errordlg('please input the starting position','parameter lost');
    return
end
if end_string==' '
    errordlg('please input the ending position','parameter lost');
    return
end
if window_string==' '
    errordlg('please input width of window','parameter lost');
    return
end
if order_string==' '
    errordlg('please input the model order','parameter lost');
    return
end

data   = load(handles.filepath);
name   = fieldnames(data);
d      = name{1};
dat    = getfield(data,d); %#ok<GFLD>
order  = str2double(order_string);
spts   = str2double(start_string);
epts   = str2double(end_string);
winlen = str2double(window_string);

% core: moving window multivariate model
mov_bi_model(dat,order,spts,epts,winlen);

% [EOF]