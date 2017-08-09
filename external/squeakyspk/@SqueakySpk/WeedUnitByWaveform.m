function WeedUnitByWaveform(SS)
% SIMPLE_GUI2 Select a data set from the pop-up menu, then
% click one of the plot-type push buttons. Clicking the button
% plots the selected data in the axes.

%  Create and then hide the GUI as it is being constructed.
f = figure('Visible','off','Position',[360,500,750,285]);

%  Construct the components.
hinstruct = uicontrol('Style','text','String','Is this average waveform Good or Bad?',...
    'Position',[50,250,500,20]);
hgood = uicontrol('Style','pushbutton','String','Good',...
    'Position',[570,220,100,40],...
    'Callback',{@goodbutton_Callback});
hbad = uicontrol('Style','pushbutton','String','Bad',...
    'Position',[570,170,100,40],...
    'Callback',{@badbutton_Callback});
hgoback = uicontrol('Style','pushbutton',...
    'String','Go Back',...
    'Position',[570,120,100,40],...
    'Callback',{@gobackbutton_Callback});
hreturn = uicontrol('Style','pushbutton',...
    'String','Return Results',...
    'Position',[570,70,100,40],...
    'Callback',{@returnbutton_Callback});

ha = axes('Units','Pixels','Position',[50,60,500,185]);
align([hgood,hbad,hgoback,hreturn],'Center','None');

% Initialize the GUI.
% Change units to normalized so components resize
% automatically.
set([f,ha,hinstruct,hgood,hbad,hreturn],...
    'Units','normalized');

%Storage for bad units
badunit = [];
goodunit = [];

% Waveform Index
waveind = 1;
current_data = SS.avgwaveform(:,waveind);
hold on
plot(current_data,'k','LineWidth',2);
plot([1 length(current_data)],[0 0],'k--')
xlabel('T (samples)')
ylabel('uV')
text(length(current_data)-10,min(current_data),['# ' num2str(waveind) ' of ' num2str(size(SS.avgwaveform,2))])
axis([1 length(current_data) min(current_data)-2 max(current_data)+2]);
% Assign the GUI a name to appear in the window title.
set(f,'Name','Supervised Unit Deletion by Average Waveform')
% Move the GUI to the center of the screen.
movegui(f,'center')
% Make the GUI visible.
set(f,'Visible','on');

%  Callbacks for simple_gui. These callbacks automatically
%  have access to component handles and initialized data
%  because they are nested at a lower level.

% Push button callbacks. Each callback plots current_data in
% the specified plot type.

    function goodbutton_Callback(source,eventdata)
        
        % Append goodunit vector
        goodunit = [goodunit waveind];
        
        if waveind+1 == size(SS.avgwaveform,2)
            returnbutton_Callback()
        else
            % Plot the next unit waveform
            cla
            waveind = waveind+1;
            current_data = SS.avgwaveform(:,waveind);
            hold on
            plot(current_data,'k','LineWidth',2);
            plot([1 length(current_data)],[0 0],'k--')
            xlabel('T (samples)')
            ylabel('uV')
            plot(current_data,'k','LineWidth',2);
            xlabel('T (samples)')
            ylabel('uV')
            text(length(current_data)-10,min(current_data),['# ' num2str(waveind) ' of ' num2str(size(SS.avgwaveform,2))])
            axis([1 length(current_data) min(current_data)-2 max(current_data)+2]);
        end
        
    end

    function badbutton_Callback(source,eventdata)
        
        % Append badunit vector
        badunit = [badunit waveind];
        
        if waveind+1 == size(SS.avgwaveform,2)
            returnbutton_Callback()
        else
            % Plot the next unit waveform
            cla
            waveind = waveind+1;
            current_data = SS.avgwaveform(:,waveind);
            hold on
            plot(current_data,'k','LineWidth',2);
            plot([1 length(current_data)],[0 0],'k--')
            xlabel('T (samples)')
            ylabel('uV')
            plot(current_data,'k','LineWidth',2);
            xlabel('T (samples)')
            ylabel('uV')
            text(length(current_data)-10,min(current_data),['# ' num2str(waveind) ' of ' num2str(size(SS.avgwaveform,2))])
            axis([1 length(current_data) min(current_data)-2 max(current_data)+2]);
        end
        
    end

    function gobackbutton_Callback(source,eventdata)
        
        if waveind > 1
            waveind = waveind - 1;
            badunit(badunit==waveind) = [];
            goodunit(badunit==waveind) = [];
            
            % Plot the old unit waveform
            cla
            current_data = SS.avgwaveform(:,waveind);
            hold on
            plot(current_data,'k','LineWidth',2);
            plot([1 length(current_data)],[0 0],'k--')
            xlabel('T (samples)')
            ylabel('uV')
            plot(current_data,'k','LineWidth',2);
            xlabel('T (samples)')
            ylabel('uV')
            text(length(current_data)-10,min(current_data),['# ' num2str(waveind) ' of ' num2str(size(SS.avgwaveform,2))])
            axis([1 length(current_data) min(current_data)-2 max(current_data)+2]);
        else
            disp('Cannot go back because you are at the first waveform!')
        end
        
    end

    function returnbutton_Callback(source,eventdata)
        
        f_end = figure('Visible','off','Position',[360,500,340,100]);
        
        %  Construct the components.
        htext = uicontrol('Style','text','String','Enter Selections into SS object?',...
            'Position',[10,70,320,20]);
        hyes = uicontrol('Style','pushbutton','String','Enter Selections',...
            'Position',[10,20,150,40],...
            'Callback',{@enterbutton_Callback});
        hno = uicontrol('Style','pushbutton','String','Ignore Selections',...
            'Position',[180,20,150,40],...
            'Callback',{@quitbutton_Callback});
        align([hgood,hbad,hgoback,hreturn],'Center','None');
        set(f_end ,'Name','Terminate Supervised Unit Deletion')
        % Move the GUI to the center of the screen.
        movegui(f_end ,'center')
        % Make the GUI visible.
        set(f_end ,'Visible','on');
        
    end


    function enterbutton_Callback(source,eventdata)
        % Display contour plot of the currently selected data.
        %  Create and then hide the GUI as it is being constructed.
        SS.RemoveUnit(badunit);
        SS.methodlog = [SS.methodlog '<WeedUnitbyWaveform>'];
        close all
        SS.PlotAvgWaveform();
        return
    end

    function quitbutton_Callback(source,eventdata)
        % Display contour plot of the currently selected data.
        %  Create and then hide the GUI as it is being constructed.
        close all
        return
    end

end