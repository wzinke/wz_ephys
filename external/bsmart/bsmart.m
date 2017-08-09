function varargout = bsmart(varargin)
% MAR M-file for mar.fig
%      MAR, by itself, creates a new MAR or raises the existing
%      singleton*.
%
%      H = MAR returns the handle to a new MAR or the handle to
%      the existing singleton*.
%
%      MAR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAR.M with the given input arguments.
%
%      MAR('Property','Value',...) creates a new MAR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mar_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mar_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Copyright (c) 2006-2007 BSMART Group
% by Richard Cui
% $Revision: 0.2$ $Date: 14-Sep-2007 11:11:11$
% SHIS UT-Houston, Houston, TX 77030, USA.
%
% Lei Xu, Hualou Liang

% Edit the above text to modify the response to help mar

% Last Modified by GUIDE v2.5 16-Sep-2007 20:35:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @mar_OpeningFcn, ...
    'gui_OutputFcn',  @mar_OutputFcn, ...
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


% --- Executes just before mar is made visible.
function mar_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mar (see VARARGIN)

% Choose default command line output for mar
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% draw bsmart face figures
axes(handles.bsmartface) % Select the proper axes
imshow('bsmartface.jpg');

% UIWAIT makes mar wait for user response (see UIRESUME)
% uiwait(handles.bsmart);


% --- Outputs from this function are returned to the command line.
function varargout = mar_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


hToolbar=uitoolbar(... % Toolbar for Open and Print buttons
    'Parent',hObject,...
    'HandleVisibility','callback');
hOpenPushtool = uipushtool(... %Open toolbar button
    'Parent',hToolbar,...
    'TooltipString','Open File',...
    'CData',iconRead(fullfile(matlabroot,...
    'toolbox/matlab/icons/opendoc.mat')),...
    'HandleVisibility','callback',...
    'ClickedCallback',@Open_menu_Callback);
hPrintPushtool = uipushtool(... % Print toolbar button
    'Parent',hToolbar,...
    'TooltipString','Print Figure',...
    'CData',iconRead(fullfile(matlabroot,...
    'toolbox/matlab/icons/printdoc.mat')),...
    'HandleVisibility','callback',...
    'ClickedCallback',@Print_menu_Callback);

% --------------------------------------------------------------------
function import_data_menu_Callback(hObject, eventdata, handles)
% hObject    handle to import_data_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
readdata;

% --------------------------------------------------------------------
function Data_write_menu_Callback(hObject, eventdata, handles)
% hObject    handle to Data_write_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
writedata;

% --------------------------------------------------------------------
function Print_menu_Callback(hObject, eventdata, handles)
% hObject    handle to Print_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Callback function run when the Print menu item is selected
printdlg();


% --------------------------------------------------------------------
function exitbsmart_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to exitbsmart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Callback function run when the Close menu item is selected

% close confirmation
% selection=...
%     questdlg(['Close ' get(handles.bsmart,'Name') '?'],...
%     ['Close ' get(handles.bsmart,'Name') '?'],...
%     'Yes','No','Yes');
% if strcmp(selection,'No')
%     return;
% end

delete(handles.bsmart);


% --------------------------------------------------------------------
function aic_Callback(hObject, eventdata, handles)
% hObject    handle to aic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
aic;

% --------------------------------------------------------------------
function FFT_menu_Callback(hObject, eventdata, handles)
% hObject    handle to FFT_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mar_fft;

% --------------------------------------------------------------------
function AMAR_menu_Callback(hObject, eventdata, handles)
% hObject    handle to AMAR_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function View_menu_Callback(hObject, eventdata, handles)
% hObject    handle to View_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function File_menu_Callback(hObject, eventdata, handles)
% hObject    handle to File_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Edit_menu_Callback(hObject, eventdata, handles)
% hObject    handle to Edit_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Tools_menu_Callback(hObject, eventdata, handles)
% hObject    handle to Tools_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Plot_menu_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Help_menu_Callback(hObject, eventdata, handles)
% hObject    handle to Help_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function whiteness_Callback(hObject, eventdata, handles)
% hObject    handle to whiteness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
whiteness;


% --------------------------------------------------------------------
function consist_Callback(hObject, eventdata, handles)
% hObject    handle to consist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
consistency_test;


% --------------------------------------------------------------------
function Lyapunov_Callback(hObject, eventdata, handles)
% hObject    handle to Lyapunov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Lyapunov;


% --------------------------------------------------------------------
function whole_Callback(hObject, eventdata, handles)
% hObject    handle to whole (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
one_window;

% --------------------------------------------------------------------
function moving_Callback(hObject, eventdata, handles)
% hObject    handle to moving (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
moving_window;

% --------------------------------------------------------------------
function analysis_Callback(hObject, eventdata, handles)
% hObject    handle to analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function spectrum_Callback(hObject, eventdata, handles)
% hObject    handle to spectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
spectrum_analysis;

% --------------------------------------------------------------------
function coherence_Callback(hObject, eventdata, handles)
% hObject    handle to coherence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
coherence;

% --------------------------------------------------------------------
function granger_causality_Callback(hObject, eventdata, handles)
% hObject    handle to granger_causality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%granger_causality;


% --------------------------------------------------------------------
function Grid_View_Callback(hObject, eventdata, handles)
% hObject    handle to Grid_View (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Topo_Plotview_Callback(hObject, eventdata, handles)
% hObject    handle to Topo_Plotview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Synoptic_Plotview_Callback(hObject, eventdata, handles)
% hObject    handle to Synoptic_Plotview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Topo_Mapview_Callback(hObject, eventdata, handles)
% hObject    handle to Topo_Mapview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function view_Callback(hObject, eventdata, handles)
% hObject    handle to view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dataview;


% --------------------------------------------------------------------
function coherence_network_Callback(hObject, eventdata, handles)
% hObject    handle to coherence_network (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
coherence_network;

% --------------------------------------------------------------------
function granger_causality_network_Callback(hObject, eventdata, handles)
% hObject    handle to granger_causality_network (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
granger_causality_network;



function whole_pairwise_Callback(hObject, eventdata, handles)
% hObject    handle to whole_pairwise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

one_window_pariwise;


% --------------------------------------------------------------------
function moving_window_pairwise_Callback(hObject, eventdata, handles)
% hObject    handle to moving_window_pairwise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
moving_window_pairwise;

% --------------------------------------------------------------------
function coherence_view_Callback(hObject, eventdata, handles)
% hObject    handle to coherence_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

coherence_view;
% --------------------------------------------------------------------
function granger_causality_view_Callback(hObject, eventdata, handles)
% hObject    handle to granger_causality_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

granger_causality_view
% --------------------------------------------------------------------
function power_view_Callback(hObject, eventdata, handles)
% hObject    handle to power_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

power_view;


% --------------------------------------------------------------------
function Preprocessing_Callback(hObject, eventdata, handles)
% hObject    handle to Preprocessing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
preprocessing;


% --------------------------------------------------------------------
function power_mul_Callback(hObject, eventdata, handles)
% hObject    handle to power_mul (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
spectrum_analysis;

% --------------------------------------------------------------------
function power_bi_Callback(hObject, eventdata, handles)
% hObject    handle to power_bi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

power_pairwise
% --------------------------------------------------------------------
function Coherence_mul_Callback(hObject, eventdata, handles)
% hObject    handle to Coherence_mul (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

coherence;
% --------------------------------------------------------------------
function Granger_mul_Callback(hObject, eventdata, handles)
% hObject    handle to Granger_mul (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Granger_bi_Callback(hObject, eventdata, handles)
% hObject    handle to Granger_bi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Coherence_bi_Callback(hObject, eventdata, handles)
% hObject    handle to Coherence_bi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pairwise_coherence;

% --------------------------------------------------------------------
function GC_bi_Callback(hObject, eventdata, handles)
% hObject    handle to GC_bi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pairwise_granger_causality;

% --------------------------------------------------------------------
function gr_Callback(hObject, eventdata, handles)
% hObject    handle to gr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function test_Callback(hObject, eventdata, handles)
% hObject    handle to test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function about_bsmart_Callback(hObject, eventdata, handles)
% hObject    handle to about_bsmart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgstr = sprintf('BSMART ver 0.5 build 102\n2006 - 2007 BSMART Goup\n(Under construction)');
msgbox(msgstr,'About BSMART','help');

% --------------------------------------------------------------------
function user_manual_Callback(hObject, eventdata, handles)
% hObject    handle to user_manual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ispc
    open('Users_guide.pdf');
else
    msgbox('Please read Users_guide.PDF (under construction)','How to use BSMART','help');
end%if


% --------------------------------------------------------------------
function fun_ref_Callback(hObject, eventdata, handles)
% hObject    handle to fun_ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ispc
    open('Function_reference.pdf');
else
    msgbox('Please read Function_reference.PDF','Function Reference','help');
end%if
