function RasterWave_Comp(SS, bound, what2show, Fs)
% RASTERWAVE_COMP Rasterwave_comp method for the SqeakySpk class.
% 
% Plots the a raster image of spike times versus unit and allows
% the user to click each point to examine the raw voltage waveform it
% represents. Spikes that will currently be cleaned in the Squeaky spike
% object appear as red dots, those that will pass appear as white dots.
%
%   RASTERWAVE_COMP(SS) displays a clickable raster-plot of
%   SS properties time versus unit (or channel),
%   e.g.    SS.time = [0.1 0.21 0.22 0.9 1.1,...,N] (NX1)
%           SS.unit = [3 1 1 1 2,...,N] (NX1)
%           SS.waveform = [[],[],...,N] (MXN).
%   The lower plot waits for the user to click a point on the raster plot
%   an then plots the corresponding raw voltage trace taken from the
%   SS.waveform matrix.
%
%   RASTERWAVE_COMP(SS, bound, showdirty, Fs) creates a raster
%   wave plot for a portion of the full recording specified in bound = [start
%   stop] measured in whatever the units SS.time are in. what2show is a
%   string argument that can take three values ['both','clean', or
%   'dirty']. This determins what type of data is shown in the raster plot.
%   Those data that will be removed after cleaning ('dirty') those data
%   that have survived cleaning ('clean') or both types in different
%   colors. Fs is the sampling frequency of the waveforms provided. Default
%   value is 25 KHz.
%
%   Created by: Jon Newman (jnewman6 at gatech dot edu)
%   Location: The Georgia Institute of Technology
%   Created on: July 30, 2009
%   Last modified: Feb 09, 2010
%
%   Licensed under the GPL: http://www.gnu.org/licenses/gpl.txt
%

% check number and type of arguments
if nargin < 4 || isempty(Fs)
    Fs = 25000; % Default sampling frequecy (Hz)
end
if nargin < 3 || isempty(what2show)
    what2show = 'both'; % Default sampling frequecy (Hz)
end
if nargin < 2 || isempty(bound)
    bound = [0 SS.time(end)];
end
if nargin < 1
    error('Need to act on SqueakSpk Object');
end

% find spike times within bound defined by user
startind = find(SS.time > bound(1),1);
if bound(2) == SS.time(end)
    endind = length(SS.time);
else
    endind = find(SS.time > bound(2),1) - 1;
end

spkinterest = SS.time(startind:endind);
cleaninterest = SS.clean(startind:endind);
waveinterest = SS.waveform(:,startind:endind);
% Has the user sorted yet?
usechan = isempty(SS.unit);
if usechan
    unitinterest = SS.channel(startind:endind);
else
    unitinterest = SS.unit(startind:endind);
end

% Set up figure and plot raster array
fh = figure();
set(fh,'color','k'); % sets the background color to black
subplot(211)
hold on

% Plot the raster for spikes in the time bound of interest
switch what2show
    case 'both'
        plot(spkinterest(cleaninterest),unitinterest(cleaninterest),'w.','MarkerSize',4);
        plot(spkinterest(~cleaninterest),unitinterest(~cleaninterest),'r.','MarkerSize',4);
    case 'clean'
        spkinterest = spkinterest(cleaninterest);
        unitinterest = unitinterest(cleaninterest);
        waveinterest = waveinterest(:,cleaninterest);
        plot(spkinterest,unitinterest,'w.','MarkerSize',4);
    case 'dirty'
        spkinterest = spkinterest(~cleaninterest);
        unitinterest = unitinterest(~cleaninterest);
        waveinterest = waveinterest(:,~cleaninterest);
        plot(spkinterest,unitinterest,'r.','MarkerSize',4);
end

axis tight
h = gca;
set(h,'color','k','XColor',[1 1 1],'YColor',[1 1 1],'ylim',[0 max(unitinterest)+1])
title('Spike Raster','fontsize',13)
set(get(h,'Title'),'Color','white')
xlabel('\textbf{time (sec)}','fontsize',13,'Interpreter','Latex')
if usechan
    ylabel('\textbf{channel}','fontsize',13,'Interpreter','Latex')
else
    ylabel('\textbf{unit}','fontsize',13,'Interpreter','Latex')
end
hold off


%Create waveform plot. This requires the user to click the raster plot
%above. It finds the spike closest to the cursor and plots the
%correspoinding voltage waveform.
lowerplot = subplot(212);
xlim([1/25 3]);
set(lowerplot,'color','k','XColor',[1 1 1],'YColor',[1 1 1])
xlabel('\textbf{time (ms)}','fontsize',13,'Interpreter','Latex')
ylabel('\textbf{V ($\mu$V)}','fontsize',13,'Interpreter','Latex')
title('Pick a point above to view voltage waveform...','fontsize',13)
set(get(lowerplot,'Title'),'Color','white')
datratio = get(h,'DataAspectRatio');
but = 1;
while but == 1
    [xi,yi] = getclick(1); %get point on current axes from mouse click
    if ~isempty(xi)
        xd = 3.9*(spkinterest - xi)/datratio(1);
        yd = (unitinterest - yi)/datratio(2);
        [temp, ind] = min(sqrt(xd.*xd + yd.*yd));
        if exist('h1','var')
            delete(h1);
            delete(h2);
        end
        subplot(211)
        hold on
        [h1 h2] = xmark([spkinterest(ind) unitinterest(ind)],range(spkinterest)/300,datratio);
        set(h1,'color','y','linewidth',1.2)
        set(h2,'color','y','linewidth',1.2)
        subplot(212)
        t = 1000/Fs.*(1:size(waveinterest,1));
        v = waveinterest(:,ind);
        plot(t,v,'y','linewidth',3)
        set(gca,'color','k','XColor',[1 1 1],'YColor',[1 1 1])
        xlabel('\textbf{time (ms)}','fontsize',13,'Interpreter','Latex')
        ylabel('\textbf{V ($\mu$V)}','fontsize',13,'Interpreter','Latex')
        title('Pick a point above to view voltage waveform...','fontsize',13)
        set(get(gca,'Title'),'Color','white')
        xlim([1000/Fs (1000/Fs)*size(SS.waveform,1)]) % x-axis in msec
    else
        return
    end
end

% 'getclick' is needed above
    function [out1,out2,out3] = getclick(arg1)
        %GINPUT Graphical input from mouse.
        %   [X,Y] = GINPUT(N) gets N points from the current axes and returns
        %   the X- and Y-coordinates in length N vectors X and Y.  The cursor
        %   can be positioned using a mouse.  Data points are entered by pressing
        %   a mouse button or any key on the keyboard except carriage return,
        %   which terminates the input before N points are entered.
        %
        %   [X,Y] = GINPUT gathers an unlimited number of points until the
        %   return key is pressed.
        %
        %   [X,Y,BUTTON] = GINPUT(N) returns a third result, BUTTON, that
        %   contains a vector of integers specifying which mouse button was
        %   used (1,2,3 from left) or ASCII numbers if a key on the keyboard
        %   was used.
        %
        %   Examples:
        %       [x,y] = ginput;
        %
        %       [x,y] = ginput(5);
        %
        %       [x, y, button] = ginput(1);
        %
        %   See also GTEXT, UIRESTORE, UISUSPEND, WAITFORBUTTONPRESS.
        
        %   Copyright 1984-2007 The MathWorks, Inc.
        %   $Revision: 5.32.4.10 $  $Date: 2007/03/03 04:45:45 $
        
        out1 = []; out2 = []; out3 = []; y = [];
        c = computer;
        if ~strcmp(c(1:2),'PC')
            tp = get(0,'TerminalProtocol');
        else
            tp = 'micro';
        end
        
        if ~strcmp(tp,'none') && ~strcmp(tp,'x') && ~strcmp(tp,'micro'),
            if nargout == 1,
                if nargin == 1,
                    out1 = trmginput(arg1);
                else
                    out1 = trmginput;
                end
            elseif nargout == 2 || nargout == 0,
                if nargin == 1,
                    [out1,out2] = trmginput(arg1);
                else
                    [out1,out2] = trmginput;
                end
                if  nargout == 0
                    out1 = [ out1 out2 ];
                end
            elseif nargout == 3,
                if nargin == 1,
                    [out1,out2,out3] = trmginput(arg1);
                else
                    [out1,out2,out3] = trmginput;
                end
            end
        else
            
            fig = gcf;
            figure(gcf);
            
            if nargin == 0
                how_many = -1;
                b = [];
            else
                how_many = arg1;
                b = [];
                if  ischar(how_many) ...
                        || size(how_many,1) ~= 1 || size(how_many,2) ~= 1 ...
                        || ~(fix(how_many) == how_many) ...
                        || how_many < 0
                    error('MATLAB:ginput:NeedPositiveInt', 'Requires a positive integer.')
                end
                if how_many == 0
                    ptr_fig = 0;
                    while(ptr_fig ~= fig)
                        ptr_fig = get(0,'PointerWindow');
                    end
                    scrn_pt = get(0,'PointerLocation');
                    loc = get(fig,'Position');
                    pt = [scrn_pt(1) - loc(1), scrn_pt(2) - loc(2)];
                    out1 = pt(1); y = pt(2);
                elseif how_many < 0
                    error('MATLAB:ginput:InvalidArgument', 'Argument must be a positive integer.')
                end
            end
            
            % Suspend figure functions
            state = uisuspend(fig);
            
            toolbar = findobj(allchild(fig),'flat','Type','uitoolbar');
            if ~isempty(toolbar)
                ptButtons = [uigettool(toolbar,'Plottools.PlottoolsOff'), ...
                    uigettool(toolbar,'Plottools.PlottoolsOn')];
                ptState = get (ptButtons,'Enable');
                set (ptButtons,'Enable','off');
            end
            
            set(fig,'pointer','cross');
            fig_units = get(fig,'units');
            char = 0;
            
            % We need to pump the event queue on unix
            % before calling WAITFORBUTTONPRESS
            drawnow
            
            while how_many ~= 0
                % Use no-side effect WAITFORBUTTONPRESS
                waserr = 0;
                try
                    keydown = wfbp;
                catch
                    waserr = 1;
                end
                if(waserr == 1)
                    if(ishandle(fig))
                        set(fig,'units',fig_units);
                        uirestore(state);
                        error('MATLAB:ginput:Interrupted', 'Interrupted');
                    else
                        return
                    end
                end
                
                ptr_fig = get(0,'CurrentFigure');
                if(ptr_fig == fig)
                    if keydown
                        char = get(fig, 'CurrentCharacter');
                        button = abs(get(fig, 'CurrentCharacter'));
                        scrn_pt = get(0, 'PointerLocation');
                        set(fig,'units','pixels')
                        loc = get(fig, 'Position');
                        % We need to compensate for an off-by-one error:
                        pt = [scrn_pt(1) - loc(1) + 1, scrn_pt(2) - loc(2) + 1];
                        set(fig,'CurrentPoint',pt);
                    else
                        button = get(fig, 'SelectionType');
                        if strcmp(button,'open')
                            button = 1;
                        elseif strcmp(button,'normal')
                            button = 1;
                        elseif strcmp(button,'extend')
                            button = 2;
                        elseif strcmp(button,'alt')
                            button = 3;
                        else
                            error('MATLAB:ginput:InvalidSelection', 'Invalid mouse selection.')
                        end
                    end
                    pt = get(gca, 'CurrentPoint');
                    
                    how_many = how_many - 1;
                    
                    if(char == 13) % & how_many ~= 0)
                        % if the return key was pressed, char will == 13,
                        % and that's our signal to break out of here whether
                        % or not we have collected all the requested data
                        % points.
                        % If this was an early breakout, don't include
                        % the <Return> key info in the return arrays.
                        % We will no longer count it if it's the last input.
                        break;
                    end
                    
                    out1 = [out1;pt(1,1)];
                    y = [y;pt(1,2)];
                    b = [b;button];
                end
            end
            
            uirestore(state);
            if ~isempty(toolbar) && ~isempty(ptButtons)
                set (ptButtons(1),'Enable',ptState{1});
                set (ptButtons(2),'Enable',ptState{2});
            end
            set(fig,'units',fig_units);
            
            if nargout > 1
                out2 = y;
                if nargout > 2
                    out3 = b;
                end
            else
                out1 = [out1 y];
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function key = wfbp
            %WFBP   Replacement for WAITFORBUTTONPRESS that has no side effects.
            
            fig = gcf;
            current_char = [];
            
            % Now wait for that buttonpress, and check for error conditions
            waserr = 0;
            try
                h=findall(fig,'type','uimenu','accel','C');   % Disabling ^C for edit menu so the only ^C is for
                set(h,'accel','');                            % interrupting the function.
                keydown = waitforbuttonpress;
                current_char = double(get(fig,'CurrentCharacter')); % Capturing the character.
                if~isempty(current_char) && (keydown == 1)           % If the character was generated by the
                    if(current_char == 3)                       % current keypress AND is ^C, set 'waserr'to 1
                        waserr = 1;                             % so that it errors out.
                    end
                end
                
                set(h,'accel','C');                                 % Set back the accelerator for edit menu.
            catch
                waserr = 1;
            end
            drawnow;
            if(waserr == 1)
                set(h,'accel','C');                                % Set back the accelerator if it errored out.
                error('MATLAB:ginput:Interrupted', 'Interrupted');
            end
            
            if nargout>0, key = keydown; end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
    end
% 'xmark' is needed above
    function [H1 H2]=xmark(center,rad,aspectR)
        scalefac = 3.9*(aspectR(2)/aspectR(1));
        X1 = center(1) + [-rad rad];
        Y1 = center(2) + [rad -rad]*scalefac;
        Y2 = center(2) + [-rad rad]*scalefac;
        H1 = line(X1,Y1);
        H2 = line(X1,Y2);
    end

end
