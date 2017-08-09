function out = ISIfilter()

global clustdata;
global clustattrib;

clustnum = inputdlg('Which cluster?','Filter non-complex spikes');
clustnum = str2num(clustnum{1});

if ~isempty(clustnum)
  
    out = false(length(clustdata.params(:,1)),1);
    index = clustattrib.clusters{clustnum}.index;
    
    
    fighandle = figure('MenuBar','none','Color',get(0,'DefaultUicontrolBackgroundColor'),'Position',[400 400 400 400], ... 
        'Tag','eventviewer','NumberTitle','off','Resize','off','Name','Choose ISI cutoff' , ... 
        'CloseRequestFcn','ISIfilter_controls(''CloseRequestFnc'',guidata(gcbo))', ...
        'WindowStyle','normal');

    
    a(1) = axes('Tag','axes1','Position',[.01 .63 .25 .25],'XTick',[],'YTick',[],'NextPlot','replacechildren','Color',[0 0 0]);
    a(2) = axes('Tag','axes2','Position',[.01 .37 .25 .25],'XTick',[],'YTick',[],'NextPlot','replacechildren','Color',[0 0 0]);
    a(3) = axes('Tag','axes3','Position',[.01 .11 .25 .25],'XTick',[],'YTick',[],'NextPlot','replacechildren','Color',[0 0 0]);
    a(4) = axes('Tag','axes4','Position',[.27 .63 .25 .25],'XTick',[],'YTick',[],'NextPlot','replacechildren','Color',[0 0 0]);
    a(5) = axes('Tag','axes5','Position',[.27 .37 .25 .25],'XTick',[],'YTick',[],'NextPlot','replacechildren','Color',[0 0 0]);
    a(6) = axes('Tag','axes6','Position',[.27 .11 .25 .25],'XTick',[],'YTick',[],'NextPlot','replacechildren','Color',[0 0 0]);

    s(1) = uicontrol('Tag','slider1','Style','slider','Units','normalized','Position',[.6 .11 .03 .77],'Callback','ISIfilter_controls(''slider1_Callback'',guidata(gcbo))','Value',1);

    t(1) = uicontrol('Tag','text1','Style','text','Units','normalized','FontUnits','pixels','FontSize',11,'Position',[.8 .81 .15 .07],'String',[num2str(1000),' ms'],'HorizontalAlignment','right');
    t(2) = uicontrol('Tag','text2','Style','text','Units','normalized','FontUnits','pixels','FontSize',11,'Position',[.67 .81 .1 .07],'String','ISI ','HorizontalAlignment','left');

    axes(a(1));
    p(1) = plot(clustdata.params(index,2),clustdata.params(index,3),'Marker','.','MarkerSize',1,'LineStyle','none','Color',[1 1 1]);
    v{1} = axis;
    hold on;
    axes(a(2));
    p(2) = plot(clustdata.params(index,2),clustdata.params(index,4),'Marker','.','MarkerSize',1,'LineStyle','none','Color',[1 1 1]);
    v{2} = axis;
    hold on;
    axes(a(3));
    p(3) = plot(clustdata.params(index,2),clustdata.params(index,5),'Marker','.','MarkerSize',1,'LineStyle','none','Color',[1 1 1]);
    v{3} = axis;
    hold on;
    axes(a(4));
    p(4) = plot(clustdata.params(index,3),clustdata.params(index,4),'Marker','.','MarkerSize',1,'LineStyle','none','Color',[1 1 1]);
    v{4} = axis;
    hold on;
    axes(a(5));
    p(5) = plot(clustdata.params(index,3),clustdata.params(index,5),'Marker','.','MarkerSize',1,'LineStyle','none','Color',[1 1 1]);
    v{5} = axis;
    hold on;
    axes(a(6));
    p(6) = plot(clustdata.params(index,4),clustdata.params(index,5),'Marker','.','MarkerSize',1,'LineStyle','none','Color',[1 1 1]);
    v{6} = axis;
    hold on;
%     %create the edit boxes
%     e(1) = uicontrol('Tag','pathedit','Style','edit','Units','normalized','HorizontalAlignment','left','BackgroundColor',[1 1 1],'Position',[.02 .85 .96 .07],'String',pwd,'Callback','ISIfilter_controls(''pathedit_Callback'',guidata(gcbo))');
%     e(2) = uicontrol('Tag','fileedit','Style','edit','Units','normalized','HorizontalAlignment','left','BackgroundColor',[1 1 1],'Position',[.15 .2 .83 .07],'String',suggestion,'Callback','ISIfilter_controls(''fileedit_Callback'',guidata(gcbo))');
%     e(2) = uicontrol('Tag','filteredit','Style','edit','Units','normalized','HorizontalAlignment','left','BackgroundColor',[1 1 1],'Position',[.15 .1 .35 .07],'String',filter,'Callback','ISIfilter_controls(''filteredit_Callback'',guidata(gcbo))');
% 
%     %create text boxes
%     t(1) = uicontrol('Tag','text1','Style','text','Units','normalized','FontUnits','pixels','FontSize',11,'Position',[.02 .2 .1 .07],'String','Filename:','HorizontalAlignment','right');
%     t(2) = uicontrol('Tag','text2','Style','text','Units','normalized','FontUnits','pixels','FontSize',11,'Position',[.02 .1 .1 .07],'String','Filter:','HorizontalAlignment','right');
% 
%     %create buttons
     b(1) = uicontrol('Tag','okbutton','Style','pushbutton','Units','normalized','Position',[.88 .03 .10 .07],'String','Ok','Callback','ISIfilter_controls(''okbutton_Callback'',guidata(gcbo))');
     b(2) = uicontrol('Tag','cancelbutton','Style','pushbutton','Units','normalized','Position',[.775 .03 .10 .07],'String','Cancel','Callback','ISIfilter_controls(''cancelbutton_Callback'',guidata(gcbo))');
     b(3) = uicontrol('Tag','signbutton','Style','pushbutton','Units','normalized','Position',[.73 .84 .08 .05],'String','>','Callback','ISIfilter_controls(''signbutton_Callback'',guidata(gcbo))','UserData',1);
%     %create the listbox
%     listcontrol(1) = uicontrol('Style','listbox','Tag','listbox1','BackgroundColor',[1 1 1],'ForegroundColor',[0 0 0], ...
%         'Units','normalized','Position',[.02 .3 .96 .5],'String',dirlist,'Max',1,'Min',1, ...
%         'Callback','ISIfilter_controls(''listbox1_Callback'',guidata(gcbo))', 'Value',[1]);


    handles = guihandles(fighandle);
    handles.index = index;
    handles.filterindex = index;
    handles.p = p;
    handles.axes = a;
    handles.axis = v;
    handles.fighandle = fighandle;
    guidata(fighandle, handles);
    
    
    for i = 1:6
        axes(handles.axes(i));
        axis(handles.axis{i});
    end



    %the gui will now force matlab to wait until one of the helper functions
    %releases it
    uiwait(fighandle);

    %when the gui is released, we output the pathname and filename, and finally
    %the figure is deleted
    handles = guidata(fighandle);
   
    delete(fighandle);

    out(handles.filterindex) = true;
    if (isempty(handles.filterindex))
        out = [];
    end



else
    out = [];
end