function varargout = one_window_pariwise(varargin)
% ONE_WINDOW_PARIWISE M-file for one_window_pariwise.fig
%      ONE_WINDOW_PARIWISE, by itself, creates a new ONE_WINDOW_PARIWISE or raises the existing
%      singleton*.
%
%      H = ONE_WINDOW_PARIWISE returns the handle to a new ONE_WINDOW_PARIWISE or the handle to
%      the existing singleton*.
%
%      ONE_WINDOW_PARIWISE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ONE_WINDOW_PARIWISE.M with the given input arguments.
%
%      ONE_WINDOW_PARIWISE('Property','Value',...) creates a new ONE_WINDOW_PARIWISE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before one_window_pariwise_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to one_window_pariwise_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright (c) 2006-2007 BSMART Group
% by Richard Cui
% $Revision: 0.3$ $Date: 14-Sep-2007 19:37:56$
% SHIS UT-Houston, Houston, TX 77030, USA.
% 
% Lei Xu, Hualou Liang


% Edit the above text to modify the response to help one_window_pariwise

% Last Modified by GUIDE v2.5 14-Sep-2007 17:55:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @one_window_pariwise_OpeningFcn, ...
                   'gui_OutputFcn',  @one_window_pariwise_OutputFcn, ...
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


% --- Executes just before one_window_pariwise is made visible.
function one_window_pariwise_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to one_window_pariwise (see VARARGIN)

% Choose default command line output for one_window_pariwise
handles.output = hObject;

% Update handles structure
handles.nchannels=0;
handles.ntrails=0;
handles.npoints=0;
handles.window=0;
handles.order=0;
handles.filepath='';
%handles.fs=0;
%handles.n=0;
%handles.freq=[];
guidata(hObject, handles);

% UIWAIT makes one_window_pariwise wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = one_window_pariwise_OutputFcn(hObject, eventdata, handles) 
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
% function pushbutton2_Callback(hObject, eventdata, handles)
% % hObject    handle to pushbutton2 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% order_string=get(handles.morder,'String');
% %fs_string=get(handles.edit3,'String');
% %n_string=get(handles.edit4,'String');
% start_string=get(handles.startpts,'String');
% window_string=get(handles.winlen,'String');
% if isequal(handles.filepath,'')
% errordlg('File not found','File Error');
% return;
% end
% data=load(handles.filepath);
% name=fieldnames(data);
% d=name{1,1};
% dat=getfield(data,d);
% si=size(dat);
%    if order_string==' '
%         errordlg('please input the model order','parameter lost');return; end
%    if fs_string==' '
%         errordlg('please input sample frequency','parameter lost');return; end
%    if start_string==' '
%         errordlg('please input the starting position','parameter lost');return; end
%    if window_string==' '
%         errordlg('please input the window','parameter lost');return; end
%    if n_string==' '
%         errordlg('please input the length of returned spectrum','parameter lost');return; end
%    handles.fs=str2num(get(handles.edit3,'String'));
%    handles.n=str2num(get(handles.edit4,'String'));
%    handles.npoints=si(1);
%    handles.ntrails=si(3);
%    handles.nchannels=si(2);
%    handles.start=str2num(get(handles.startpts,'String'));
%    handles.window=str2num(get(handles.winlen,'String'));
%    handles.order=str2num(get(handles.morder,'String'));
%    guidata(hObject,handles);
%    n=handles.n;
%    fs=handles.fs;
%    if handles.order>20
%        errordlg('please input the model order less than 20','parameter lost');return; end
%    if handles.start>handles.npoints
%        errordlg('please input the correct starting position','parameter lost');return; end
%    if handles.window>handles.npoints-handles.start+1
%        errordlg('please input the correct window length','parameter lost');return; end
%    channel2=handles.nchannels;
%    channel=2;
%    save channel channel -ascii;
%    
%    trail=handles.ntrails;
%    save trail trail -ascii;
%    
%    points=handles.window;
%    save points points -ascii;
%    
%    order=handles.order;
%    save order order -ascii;
%    dat=dat(handles.start:handles.start+handles.window-1,:,:);
%    pairco=zeros(1,n,channel2*(channel2-1)/2);
%    spect=zeros(1,n,channel2);
%    k=0;l=1;
%    for i=1:(channel2-1)
%        for j=(i+1):channel2
%            dat1=dat(:,i,:);
%            dat2=dat(:,j,:);
%            dat3=cat(2,dat1,dat2);
%            writedat('dataset.bin',dat3);
%            eval(['unix ' '(''' 'opssfull ' 'dataset.bin ' ' A ' 'Ve ' 'AIC' ''')']);
%            aredat=load('A');
%            arndat=load('Ve');
%            [csd4d, freq] = MAR_csd4dmem(aredat,arndat,n,fs);
%            [paircoh, partcoh, mulcoh, autospect] = MAR_coh3D(csd4d, fs);
%            k=k+1;
%            pairco(1,:,k)=paircoh(1,:,1);
% 
%            if i==1
%                l=l+1;
%            spect(1,:,i)=autospect(1,:,1);
%            spect(1,:,l)=autospect(1,:,2);
%     
%            end
%            !DEL A
%            !DEL Ve
%            !DEL AIC
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
% function pushbutton4_Callback(hObject, eventdata, handles)
% % hObject    handle to pushbutton4 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% global Fx2y Fy2x channel;
% order_string=get(handles.morder,'String');
% fs_string=get(handles.edit3,'String');
% start_string=get(handles.startpts,'String');
% window_string=get(handles.winlen,'String');
% freq_string=get(handles.edit5,'String');
% 
% if isequal(handles.filepath,'')
% errordlg('File not found','File Error');
% return;
% end
% data=load(handles.filepath);
% name=fieldnames(data);
% d=name{1,1};
% dat=getfield(data,d);
% si=size(dat);
%    if order_string==' '
%         errordlg('please input the model order','parameter lost');return; end
%    if fs_string==' '
%         errordlg('please input sample frequency','parameter lost');return; end
%    if start_string==' '
%         errordlg('please input the starting position','parameter lost');return; end
%    if window_string==' '
%         errordlg('please input the window','parameter lost');return; end
%    if freq_string==' '
%         errordlg('please input a vector of frequencies of interest, usually freq=0:fs/2','parameter lost');return; end
%    handles.fs=str2num(get(handles.edit3,'String'));
%    handles.npoints=si(1);
%    handles.ntrails=si(3);
%    handles.nchannels=si(2);
%    handles.freq=eval(get(handles.edit5,'String'));
%    handles.order=str2num(get(handles.morder,'String'));
%    handles.start=str2num(get(handles.startpts,'String'));
%    handles.window=str2num(get(handles.winlen,'String'));
%    guidata(hObject,handles);
%    fs=handles.fs;
%    if handles.order>20
%        errordlg('please input the model order less than 20','parameter lost');return; end
%    if handles.start>handles.npoints
%        errordlg('please input the correct starting position','parameter lost');return; end
%    if handles.window>handles.npoints-handles.start+1
%        errordlg('please input the correct window length','parameter lost');return; end
%    channel2=handles.nchannels;
%    channel=channel2;
%    trail=handles.ntrails;
%    points=handles.window;
%    order=handles.order;
%    fs=handles.fs;
%    freq=handles.freq;
%    a=zeros(2,trail*points);
%    dat=dat(handles.start:handles.start+handles.window-1,:,:);
%    Fxy=[];
%    Fyx=[];
%    for i=1:(channel2-1)
%        for j=(i+1):channel2
%            for k=1:trail
%                a(1,(k-1)*points+1:k*points)=dat(:,i,k);
%                a(2,(k-1)*points+1:k*points)=dat(:,j,k);
%            end
%            [pp,cohe,Fx2y,Fy2x]=pwcausal(a,trail,points,order,fs,freq);
%            Fxy=[Fxy;Fx2y];
%            Fyx=[Fyx;Fy2x];
%        end
%    end
%    Fx2y=Fxy;
%    Fy2x=Fyx;
%    uisave Fx2y Fx2y;
%    uisave Fy2x Fy2x;
%  granger_causality_analysis;


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





function startpts_Callback(hObject, eventdata, handles)
% hObject    handle to startpts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startpts as text
%        str2double(get(hObject,'String')) returns contents of startpts as a double


% --- Executes during object creation, after setting all properties.
function startpts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startpts (see GCBO)
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


% --- Executes on button press in create_model.
function create_model_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to create_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(handles.filepath,'')
    errordlg('File not found','File Error');
    return
end

order_string  = get(handles.morder,'String');
start_string  = get(handles.startpts,'String');
window_string = get(handles.winlen,'String');
if order_string == ' '
    errordlg('please input the model order','parameter lost');
    return
end
if start_string == ' '
    errordlg('please input the starting position','parameter lost');
    return
end
if window_string == ' '
    errordlg('please input the window','parameter lost');
    return
end

data = load(handles.filepath);
name = fieldnames(data);
d = name{1};
dat    = getfield(data,d); %#ok<GFLD>
order  = str2double(order_string);
spts   = str2double(start_string);
winlen = str2double(window_string);

% core: fixed window bivariate model
h = msgbox('Please waiting...','Fixed window bivariate model','help');
one_bi_model(dat,order,spts,winlen);
close(h);

% [EOF]    
