function varargout = grid_view(varargin)
% GRID_VIEW M-file for grid_view.fig
%      GRID_VIEW, by itself, creates a new GRID_VIEW or raises the existing
%      singleton*.
%
%      H = GRID_VIEW returns the handle to a new GRID_VIEW or the handle to
%      the existing singleton*.
%
%      GRID_VIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GRID_VIEW.M with the given input arguments.
%
%      GRID_VIEW('Property','Value',...) creates a new GRID_VIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before grid_view_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to grid_view_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright (c) 2006-2007 BSMART group.
% by Richard Cui
% $Revision: 0.2$ $Date: 12-Sep-2007 16:41:40$
% SHIS UT-Houston, Houston, TX 77030, USA.
%
% Lei Xu, Hualou Liang

% Edit the above text to modify the response to help grid_view

% Last Modified by GUIDE v2.5 15-Sep-2007 23:10:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @grid_view_OpeningFcn, ...
                   'gui_OutputFcn',  @grid_view_OutputFcn, ...
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


% --- Executes just before grid_view is made visible.
function grid_view_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to grid_view (see VARARGIN)

% Choose default command line output for grid_view
handles.output = hObject;

% Update handles structure
data = varargin{1};
handles.data = data;
guidata(hObject, handles);

% UIWAIT makes grid_view wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = grid_view_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.figure1);

% --- Executes on button press in plot.
function plot_Callback(hObject, eventdata, handles)
% hObject    handle to plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% data
dat = handles.data;
npts = size(dat,1);
nchn = size(dat,2);
ntrl = size(dat,3);

% specified trial
trial_string = get(handles.trial,'String');
if trial_string==' '
       errordlg('please input the trail you interested','parameter lost');
       return
end
trial = str2double(trial_string);
if (trial > ntrl || trial <= 0)
    errordlg('please input correct trail number','parameter lost');
    return
end

% specified channels range
chstart_string = get(handles.chstart,'String');
if chstart_string==' '
       errordlg('please input start channel','parameter lost');
       return
end
chstart = str2double(chstart_string);

chend_string   = get(handles.chend,'String');
if chend_string==' '
    errordlg('please input end channel','parameter lost');
    return
end
chend = str2double(chend_string);

if (chstart > nchn || chend <= 0 || chend > nchn || chend <= 0 || chstart > chend )
    errordlg('please input correct channel number','parameter lost');
    return
end

% specified time range
spts_string  = get(handles.timestart,'String');
if spts_string==' '
    errordlg('please input started point','parameter lost');
    return
end
spts = str2double(spts_string);

epts_string  = get(handles.timeend,'String');
if epts_string==' '
    errordlg('please input ended point','parameter lost');
    return
end
epts = str2double(epts_string);

if (spts > npts || spts <= 0 || epts > npts || epts <= 0 || spts > epts )
    errordlg('please input correct point start and end location','parameter lost');
    return
end

% core: draw grid view
gridview(dat,trial,chstart,chend,spts,epts);

% 
% 
% guidata(hObject, handles);
% figure('Name','Grid View','NumberTitle','off')
%     x=data(po1:po2,ch1:ch2,handles.trail);
%     imagesc([po1:po2],[ch1:ch2],x');
%     h=gca;
%     set(h,'XLim',[po1 po2]);
%     set(h,'XTick',[po1:po2]);
%     set(h,'YTick',[ch1:ch2]);
%     xlabel('TimePoints');
%     ylabel('Channel');
% colorbar;
% figure('Name','Grid View 2','NumberTitle','off');
% for i=ch1:1:ch2 
%     subplot(ch2-ch1+1,1,i-ch1+1);
%     x=data(po1:po2,i,handles.trail);
%     plot([po1:po2],x);
%     h=gca;
%     set(h,'XLim',[po1 po2]);
%     set(h,'XTick',[po1:po2]);
%     ylabel(num2str(i));
%     if i==ch2
%         xlabel('TimePoints');
%     end
% end
%     

function chstart_Callback(hObject, eventdata, handles)
% hObject    handle to chstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of chstart as text
%        str2double(get(hObject,'String')) returns contents of chstart as a double


% --- Executes during object creation, after setting all properties.
function chstart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function chend_Callback(hObject, eventdata, handles)
% hObject    handle to chend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of chend as text
%        str2double(get(hObject,'String')) returns contents of chend as a double


% --- Executes during object creation, after setting all properties.
function chend_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function trial_Callback(hObject, eventdata, handles)
% hObject    handle to trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trial as text
%        str2double(get(hObject,'String')) returns contents of trial as a double


% --- Executes during object creation, after setting all properties.
function trial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function timeend_Callback(hObject, eventdata, handles)
% hObject    handle to timeend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of timeend as text
%        str2double(get(hObject,'String')) returns contents of timeend as a double


% --- Executes during object creation, after setting all properties.
function timeend_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timeend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function timestart_Callback(hObject, eventdata, handles)
% hObject    handle to timestart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of timestart as text
%        str2double(get(hObject,'String')) returns contents of timestart as a double


% --- Executes during object creation, after setting all properties.
function timestart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timestart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


