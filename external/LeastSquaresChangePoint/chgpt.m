function varargout = chgpt(varargin)
% CHGPT M-file for chgpt.fig
%      CHGPT, by itself, creates a new CHGPT or raises the existing
%      singleton*.
%
%      H = CHGPT returns the handle to a new CHGPT or the handle to
%      the existing singleton*.
%
%      CHGPT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHGPT.M with the given input arguments.
%
%      CHGPT('Property','Value',...) creates a new CHGPT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before chgpt_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to chgpt_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help chgpt

% Last Modified by GUIDE v2.5 25-Sep-2009 10:57:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @chgpt_OpeningFcn, ...
                   'gui_OutputFcn',  @chgpt_OutputFcn, ...
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


% --- Executes just before chgpt is made visible.
function chgpt_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output skip_checkbox, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to chgpt (see VARARGIN)

% Choose default command line output for chgpt
handles.output = hObject;

%for radio button group
set(handles.weight_buttons,'SelectionChangeFcn',@weight_buttons_SelectionChangeFcn);
set(handles.regressor_buttons,'SelectionChangeFcn',@regressor_buttons_SelectionChangeFcn);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes chgpt wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = chgpt_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output skip_checkbox (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function filename_input_Callback(hObject, eventdata, handles)
% hObject    handle to filename_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filename_input as text
%        str2double(get(hObject,'String')) returns contents of filename_input as a double
input=get(hObject, 'String');
if(~exist(input, 'file'))
    set(hObject,'String','file.txt')
    set(handles.status, 'String', 'File Does Not Exist');
end

% --- Executes during object creation, after setting all properties.
function filename_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filename_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kmax_input_Callback(hObject, eventdata, handles)
% hObject    handle to kmax_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kmax_input as text
%        str2double(get(hObject,'String')) returns contents of kmax_input as a double

%store the contents of kmax_input as a string. if the string
%is not a number then input will be empty
input = str2num(get(hObject,'String'));

%checks to see if input is empty. if so, default input1_editText to zero
if (isempty(input))
     set(hObject,'String','10')
     set(handles.status, 'String', 'Invalid Entry for Maximal Number of Change Points');
end
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function kmax_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kmax_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_chgpts_input_Callback(hObject, eventdata, handles)
% hObject    handle to num_chgpts_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_chgpts_input as text
%        str2double(get(hObject,'String')) returns contents of num_chgpts_input as a double

%store the contents of input1_editText as a string. if the string
%is not a number then input will be empty
input = str2num(get(hObject,'String'));

%checks to see if input is empty. if so, default input1_editText to zero
if (isempty(input))
     set(hObject,'String','5')
     set(handles.status, 'String', 'Invalid Entry for Number of Change Points');
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function num_chgpts_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_chgpts_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
goodtogo=1;     %have all variables been entered properly?

%retrieve GUI data, i.e. the handles structure and check for valid entries
arg=get(handles.skip_checkbox,'Value');
if (arg)
    skip='s';
else skip ='n';
end
filename=(get(handles.filename_input, 'String'));
if(~exist(filename, 'file') && skip=='n')
    set(handles.filename_input,'String','file.txt')
    set(handles.status, 'String', 'File Does Not Exist');
    goodtogo=0;
end
num_chgpts = str2num(get(handles.num_chgpts_input,'String'));
if (num_chgpts==0)
    goodtogo=0;
    set(handles.status, 'String', 'Number of Change Points for Output Must Be Larger than Zero');
end
kmax = str2num(get(handles.kmax_input,'String'));
if (kmax==0 || kmax<num_chgpts)
    goodtogo=0;
    set(handles.status, 'String', 'Maximum Number of Change Points Must Be Larger than Number of Change Points in Output');
end

pause on;
if (goodtogo)
    c = 'Running';
    % need to convert the answer back into String type to display it
    set(handles.status,'String',c);
    guidata(hObject, handles);
    pause(0.01);
    changepoint(handles);
end
guidata(hObject, handles);

% --- Executes on button press in detrend_checkbox.
function detrend_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to detrend_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of detrend_checkbox

% --- Executes on button press in skip_checkbox.
function skip_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to skip_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of skip_checkbox



function dist_input_Callback(hObject, eventdata, handles)
% hObject    handle to dist_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dist_input as text
%        str2double(get(hObject,'String')) returns contents of dist_input as a double

%store the contents of dist_input as a string. if the string
%is not a number then input will be empty
input = str2num(get(hObject,'String'));

%checks to see if input is empty. if so, default input1_editText to zero
if (isempty(input))
     set(hObject,'String','0')
     set(handles.status, 'String', 'Invalid Entry for Distance Between Change Points');
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function dist_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dist_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function weight_buttons_SelectionChangeFcn(hObject, eventdata)
 %Unnecessary??
%retrieve GUI data, i.e. the handles structure
handles = guidata(hObject); 
 
%updates the handles structure
guidata(hObject, handles);


% --- Executes on button press in browse_pushbutton.
function browse_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to browse_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname]=uigetfile('.txt');
set(handles.filename_input, 'String', [pathname filename]);
%set(handles.fname_text, 'String', filename);

function w1_Callback(hObject, eventdata, handles)
% hObject    handle to w1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of w1 as text
%        str2double(get(hObject,'String')) returns contents of w1 as a double
input = str2num(get(hObject,'String'));

%checks to see if input is empty. if so, default w1 to blank
if (isempty(input))
     set(hObject,'String','')
     set(handles.status, 'String', 'Invalid Entry for Regressor - Periodicity Must Be Larger than 0');
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function w1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to w1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function w4_Callback(hObject, eventdata, handles)
% hObject    handle to w4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of w4 as text
%        str2double(get(hObject,'String')) returns contents of w4 as a double
input = str2num(get(hObject,'String'));

%checks to see if input is empty. if so, default w1 to blank
if (isempty(input))
     set(hObject,'String','')
     set(handles.status, 'String', 'Invalid Entry for Regressor - Periodicity Must Be Larger than 0');
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function w4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to w4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function w2_Callback(hObject, eventdata, handles)
% hObject    handle to w2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of w2 as text
%        str2double(get(hObject,'String')) returns contents of w2 as a double
input = str2num(get(hObject,'String'));

%checks to see if input is empty. if so, default w1 to blank
if (isempty(input))
     set(hObject,'String','')
     set(handles.status, 'String', 'Invalid Entry for Regressor - Periodicity Must Be Larger than 0');
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function w2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to w2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function w5_Callback(hObject, eventdata, handles)
% hObject    handle to w5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of w5 as text
%        str2double(get(hObject,'String')) returns contents of w5 as a double
input = str2num(get(hObject,'String'));

%checks to see if input is empty. if so, default w1 to blank
if (isempty(input))
     set(hObject,'String','')
     set(handles.status, 'String', 'Invalid Entry for Regressor - Periodicity Must Be Larger than 0');
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function w5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to w5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function w3_Callback(hObject, eventdata, handles)
% hObject    handle to w3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of w3 as text
%        str2double(get(hObject,'String')) returns contents of w3 as a double
input = str2num(get(hObject,'String'));

%checks to see if input is empty. if so, default w1 to blank
if (isempty(input))
     set(hObject,'String','')
     set(handles.status, 'String', 'Invalid Entry for Regressor - Periodicity Must Be Larger than 0');
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function w3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to w3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function regressor_buttons_SelectionChangeFcn(hObject, eventdata)
 
%retrieve GUI data, i.e. the handles structure
handles = guidata(hObject); 
 
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'file_radiobutton'
      %execute this code when file_radiobutton is selected
      set(handles.w1,'String','');
      set(handles.w2,'String','');
      set(handles.w3,'String','');
      set(handles.w4,'String','');
      set(handles.w5,'String','');
    set(handles.poly_order, 'String', '');
    
    case 'create_radiobutton'
      %execute this code when create_radiobutton is selected
      set(handles.poly_order, 'String', '');
      
    case 'poly_radiobutton'
        set(handles.poly_order, 'String', '0');
        set(handles.w1,'String','');
      set(handles.w2,'String','');
      set(handles.w3,'String','');
      set(handles.w4,'String','');
      set(handles.w5,'String','');
    
    otherwise
       % Code for when there is no match.
 
end
%updates the handles structure
guidata(hObject, handles);



function poly_order_Callback(hObject, eventdata, handles)
% hObject    handle to poly_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of poly_order as text
%        str2double(get(hObject,'String')) returns contents of poly_order as a double
input = str2num(get(hObject,'String'));

%checks to see if input is empty. if so, default w1 to blank
if (isempty(input)||input<0 || input>3)
     set(hObject,'String','0')
     set(handles.status, 'String', 'Invalid Entry for Order of Polynomial - Must be Between 0 and 3');
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function poly_order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to poly_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


