function varargout = AIC_order_view(varargin)
% AIC_ORDER_VIEW M-file for AIC_order_view.fig
%      AIC_ORDER_VIEW, by itself, creates a new AIC_ORDER_VIEW or raises the existing
%      singleton*.
%
%      H = AIC_ORDER_VIEW returns the handle to a new AIC_ORDER_VIEW or the handle to
%      the existing singleton*.
%
%      AIC_ORDER_VIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AIC_ORDER_VIEW.M with the given input arguments.
%
%      AIC_ORDER_VIEW('Property','Value',...) creates a new AIC_ORDER_VIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AIC_order_view_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AIC_order_view_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AIC_order_view

% Last Modified by GUIDE v2.5 14-Sep-2007 17:36:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AIC_order_view_OpeningFcn, ...
                   'gui_OutputFcn',  @AIC_order_view_OutputFcn, ...
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


% --- Executes just before AIC_order_view is made visible.
function AIC_order_view_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AIC_order_view (see VARARGIN)

% Choose default command line output for AIC_order_view
handles.output = hObject;

points = varargin{1};
winlen = varargin{2};
aic    = varargin{3};
handles.points = points;
handles.winlen = winlen;
handles.AIC = aic;

% Update handles structure
handles.window = 0;
guidata(hObject, handles);

% UIWAIT makes AIC_order_view wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AIC_order_view_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in aic_plot.
function aic_plot_Callback(hObject, eventdata, handles)
% hObject    handle to aic_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

window_string = get(handles.winpos,'String');
if window_string==' '
    errordlg('please input a window number','parameter lost');
    return
end

wind = str2double(window_string);
winlen = handles.winlen;
points = handles.points;
if (wind > (points-winlen+1) || wind <= 0)
    errordlg('please input correct window number','parameter lost');
    return
end

aic = handles.AIC;
dat = aic(wind,:);
data=dat(~isnan(dat));

% draw the figure
figure('Name','AIC Model Order Estimation','NumberTitle','off')
plot(data,'-s');
h = gca;
xlabel(h,'Model Order');
ylabel(h,'AIC measure');


% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.figure1);



function winpos_Callback(hObject, eventdata, handles)
% hObject    handle to winpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of winpos as text
%        str2double(get(hObject,'String')) returns contents of winpos as a double


% --- Executes during object creation, after setting all properties.
function winpos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to winpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


