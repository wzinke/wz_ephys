function varargout = consistency_test(varargin)
%CONSISTENCY_TEST M-file for consistency_test.fig
%      CONSISTENCY_TEST, by itself, creates a new CONSISTENCY_TEST or raises the existing
%      singleton*.
%
%      H = CONSISTENCY_TEST returns the handle to a new CONSISTENCY_TEST or the handle to
%      the existing singleton*.
%
%      CONSISTENCY_TEST('Property','Value',...) creates a new CONSISTENCY_TEST using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to consistency_test_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      CONSISTENCY_TEST('CALLBACK') and CONSISTENCY_TEST('CALLBACK',hObject,...) call the
%      local function named CALLBACK in CONSISTENCY_TEST.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright (c) 2006-2007 BSMART Group
% by Richard Cui
% $Revision: 0.2$ $Date: 15-Sep-2007 09:04:13$
% SHIS UT-Houston, Houston, TX 77030, USA.
% 
% Lei Xu, Hualou Liang

% Edit the above text to modify the response to help consistency_test

% Last Modified by GUIDE v2.5 15-Sep-2007 09:02:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @consistency_test_OpeningFcn, ...
    'gui_OutputFcn',  @consistency_test_OutputFcn, ...
    'gui_LayoutFcn',  [], ...
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


% --- Executes just before consistency_test is made visible.
    function consistency_test_OpeningFcn(hObject, eventdata, handles, varargin)
        % This function has no output args, see OutputFcn.
        % hObject    handle to figure
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)
        % varargin   unrecognized PropertyName/PropertyValue pairs from the
        %            command line (see VARARGIN)

        % Choose default command line output for consistency_test
        handles.output = hObject;

        % Update handles structure
        handles.nchannels=0;
        handles.ntrails=0;
        handles.npoints=0;
        handles.window=0;
        handles.order=0;
        handles.filepath='';
        handles.arcfilepath='';
        handles.arnfilepath='';
        guidata(hObject, handles);

        % UIWAIT makes consistency_test wait for user response (see UIRESUME)
        % uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
        function varargout = consistency_test_OutputFcn(hObject, eventdata, handles)
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


% --- Executes on button press in choose_nvm.
function choose_nvm_Callback(hObject, eventdata, handles)
% hObject    handle to choose_nvm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
{'*.*',  'All Files (*.*)'; ...
'*.m;*.fig;*.mat;*.mdl','MATLAB Files (*.m,*.fig,*.mat,*.mdl)';
   '*.m',  'M-files (*.m)'; ...
   '*.fig','Figures (*.fig)'; ...
   '*.mat','MAT-files (*.mat)'; ...
   '*.mdl','Models (*.mdl)'; ...
   '*.dm6','steven file'},...
   'Pick a file');
if isequal(filename,0) | isequal(pathname,0)
else
    handles.arnfilepath=fullfile(pathname,filename);
end
    guidata(hObject,handles);
    set(handles.fnvm_show,'String',handles.arnfilepath);


% --- Executes on button press in choose_arc.
function choose_arc_Callback(hObject, eventdata, handles)
% hObject    handle to choose_arc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
{'*.*',  'All Files (*.*)'; ...
'*.m;*.fig;*.mat;*.mdl','MATLAB Files (*.m,*.fig,*.mat,*.mdl)';...
   '*.m',  'M-files (*.m)'; ...
   '*.fig','Figures (*.fig)'; ...
   '*.mat','MAT-files (*.mat)'; ...
   '*.mdl','Models (*.mdl)'; ...
   '*.dm6','steven file'},...
   'Pick a file');
if isequal(filename,0) | isequal(pathname,0)
else
handles.arcfilepath=fullfile(pathname,filename);
end
guidata(hObject,handles);
set(handles.farc_show,'String',handles.arcfilepath);



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
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function farc_show_Callback(hObject, eventdata, handles)
% hObject    handle to farc_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of farc_show as text
%        str2double(get(hObject,'String')) returns contents of farc_show as a double


% --- Executes during object creation, after setting all properties.
function farc_show_CreateFcn(hObject, eventdata, handles)
% hObject    handle to farc_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function fnvm_show_Callback(hObject, eventdata, handles)
% hObject    handle to fnvm_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fnvm_show as text
%        str2double(get(hObject,'String')) returns contents of fnvm_show as a double


% --- Executes during object creation, after setting all properties.
function fnvm_show_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fnvm_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in close.
    function close_Callback(hObject, eventdata, handles)
        % hObject    handle to close (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)
        delete(handles.consistency_figure);

        % --- Executes on button press in con_test.
        function con_test_Callback(hObject, eventdata, handles)
            % hObject    handle to con_test (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            if isequal(handles.filepath,'')
                errordlg('Data File not found','File Error');
                return
            end
            if isequal(handles.arcfilepath,'')
                errordlg('Coffecient File not found','File Error');
                return
            end
            if isequal(handles.arnfilepath,'')
                errordlg('Noise File not found','File Error');
                return
            end
           
            % data
            data  = load(handles.filepath);
            name  = fieldnames(data);
            d     = name{1};
            dat   = getfield(data,d); %#ok<GFLD>
            % AR coefficients
            arc   = load(handles.arcfilepath);
            name  = fieldnames(arc);
            d     = name{1};
            A     = getfield(arc,d); %#ok<GFLD>
            % make it compatable with Matlab R14
            if size(A,2) == 1
                A = A.';    % if A is a column, change it to row
            end%if
            % Noise variance matrix
            nvm   = load(handles.arnfilepath);
            name  = fieldnames(nvm);
            d     = name{1};
            Ve   = getfield(nvm,d); %#ok<GFLD>
            % make it compatable with Matlab R14
            if size(Ve,2) == 1
                Ve = Ve.';    % if Ve is a column, change it to row
            end%if

            % core: consistency test
            h = msgbox('Please wait...','Consistency Text','help');
            ratio = consistencytest(dat,A,Ve);
            close(h);
            
            % save results
            suggestname = 'Per_con';
            [filename,pathname] = uiputfile( ...
                {'*.mat','MAT-files (*.mat)';...
                '*.*',  'All Files (*.*)'},...
                'Save as',...
                suggestname);

            if isequal(filename,0) || isequal(pathname,0)
            else
                filepath = fullfile(pathname,filename);
                save(filepath,'ratio');
            end%if

            % plot results
            figure('Name','Consistency Test','NumberTitle','off')
            plot(ratio,'-o');
            h = gca;
            xlabel(h,'Time')
            ylabel(h,'Percent Consistency');
            set(h,'YLim',[0 100]);


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
                if isequal(filename,0) || isequal(pathname,0)
                else
                    handles.filepath=fullfile(pathname,filename);
                end
                guidata(hObject,handles);
                set(handles.fname_show,'String',handles.filepath);

