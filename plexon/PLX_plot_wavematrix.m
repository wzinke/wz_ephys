function [fhndl, mplots] = PLX_plot_wavematrix(plxfile, plot_unsorted, matsz, spkchan, do_mean, flnm)
% Plot the median spike waveforms for all channels.
%  
% DESCRIPTION 
%   This functions creates a simple figure matrix that shows the median spike
%   waveform and the 25/75-percentiles of all spike waveforms stored in a plexon
%   file.
%
%   This routine needs the readPLXFileC function provided by B. Kraus on Matlab
%   File Exchange:
%   http://www.mathworks.com/matlabcentral/fileexchange/42160-readplxfilec
%
%
% SYNTAX 
% 
%   [fhndl] = PLX_plot_wavematrix(plxfile, plot_unsorted)
%
%   Input:
%
%       plxfile         name of a plexon data file
%
%       plot_unsorted   also plot average waveform for spikes not assigned to a unit.
%
%       matsz           specify number of channels to represent. If not data is
%                       availabe for a current channel the plot will be empty. 
%                       If not specified, the matrix size could vary depending 
%                       on the numbe rof channels with data
%       
%       spkchan         plot only channels specified here
%
%       do_mean         plot mean and standard deviation insteat of median and quartiles
%       
%       flnm            filename to safe plot as image
%
%
%   Output:
%
%       fhndl           handle for the figure
%
%       mplots          vector with handles to all subplots
%
% ......................................................................... 
% wolf zinke, wolfzinke@gmail.com 
%
% $Created : 11-Feb-2015 by wolf zinke

% ____________________________________________________________________________ %
%% Set plot parameters

% define default colors, assume a fixed maximum number of spikes per channel
main_col = [0.6000, 0.6000, 0.6000; ...
                 0,      0, 1.0000; ...
            1.0000,      0,      0; ...
                 0, 1.0000,      0; ...
            1.0000, 0.1034, 0.7241; ...
            1.0000, 0.8276,      0];
ptch_col = 1 - ((1-main_col) .* 0.25);

% ____________________________________________________________________________ %
%% set default arguments
if(~exist('plot_unsorted','var') || isempty(plot_unsorted))
    plot_unsorted = 0;
end
if(~exist('do_mean','var') || isempty(do_mean))
    do_mean = 0;
end

if(~exist('flnm','var'))
    flnm = [];
end

% ____________________________________________________________________________ %
%% Check if dependencies are fullfilled
if(~exist('readPLXFileC','file'))
    warning('readPLXFileC not found, please select now');
    [~,PathName] = uigetfile({'*.mex*'},'readPLXFileC mex file');
    addpath(PathName);
end

% if not specified use GUI to get the file
if(~exist('plxfile','var') || isempty(plxfile))
    [FileName,PathName] = uigetfile({'*.plx;*.PLX'},'Load plexon file');
    plxfile = fullfile(PathName,FileName);
end

% ____________________________________________________________________________ %
%% load plexon file

disp(['Reading plexon file: ',plxfile]); tm = cputime;
plx = readPLXFileC(plxfile,'waves','noevents', 'nocontinuous');
fprintf('... done after %.1f s\n\n', cputime-tm);

[~, plxname] = fileparts(plxfile);

if(~exist('spkchan','var') || isempty(spkchan))
    unitcnt = sum(plx.SpikeTimestampCounts > 0);
    spkchan = find(unitcnt > 0);
end

numUnits = length(spkchan);

if(~exist('matsz','var') || isempty(matsz))
    matsz = numUnits;
end

if(matsz == numUnits)
    check_channel = 0; % no need to check if the channel is valid, because only valid channels will be plotted
else
    check_channel = 1; 
end

% there must be a much nicer way to code this one...
if(matsz <=2)
    num_rows = 1;
    num_cols = 2;
elseif(matsz <=4)
    num_rows = 2;
    num_cols = 2;
elseif(matsz <=6)
    num_rows = 2;
    num_cols = 3;
elseif(matsz <=9)
    num_rows = 3;
    num_cols = 3;
elseif(matsz <=6)
    num_rows = 2;
    num_cols = 3;
elseif(matsz <=12)
    num_rows = 3;
    num_cols = 4;
elseif(matsz <=16)
    num_rows = 4;
    num_cols = 4;
elseif(matsz <=20)
    num_rows = 4;
    num_cols = 5;
elseif(matsz <=24)
    num_rows = 4;
    num_cols = 6;
end   
    
figsz  = [0 0 num_cols*200 num_rows*150];
Ystart = linspace(97,3,num_rows+1)./100;
Ystart(1)   = [];
Xstart = linspace(3,97,num_cols+1)./100;
Xstart(end) = [];

xwd = min(diff(Xstart))-0.004;
ywd = min(abs(diff(Ystart)))-0.006;

% ____________________________________________________________________________ %
%% loop over channels and create a plot
chancnt = 0; % currently this assumes that channel count starts at 1. Needs to be corrected to be more flexible!

fhndl  = figure('Name', plxname, 'Position', figsz, 'Renderer', 'Painters');
mplots = []; chnlst = [];
minY   = inf; maxY = -inf;

for(x=1:num_cols)
    cX = Xstart(x);
    for(y=1:num_rows)
        cY = Ystart(y);
        
        chancnt = chancnt+1;
        
        if(chancnt > matsz)
            break
        end
        
        chndl = subplot('Position',[cX cY xwd ywd]);
        hold on;
%        hline(0);
        plot(xlim,[0,0]);
        
        if(x>1)
            set(gca,'YTickLabel',[]);
        end
    
        if(check_channel == 1)
             curr_chan = find(spkchan == chancnt);
        else
            curr_chan = spkchan(chancnt);
        end
                
        % check if current channel is in the list
        if(~isempty(curr_chan))
            
            unitlst = unique(plx.SpikeChannels(curr_chan).Units);
            
            if(plot_unsorted == 0) % remove unsorted spike wave forms
                unitlst(unitlst==0) = [];
            end
            
            if(~isempty(unitlst))
                mplots = [mplots, chndl];
                chnlst = [chnlst, curr_chan];
                
                unitwaves = double(plx.SpikeChannels(curr_chan).Waves)'; % get waveform matrix
                
                for(u=1:length(unitlst))
                    cunit = unitlst(u);
                    spkpos = plx.SpikeChannels(curr_chan).Units == cunit;
                    cwaves = unitwaves(spkpos,:)./1000;
                    
                    if(do_mean == 0)
                        spikewave = prctile(cwaves, 50);
                        h_val     = prctile(cwaves, 75);
                        l_val     = prctile(cwaves, 25);
                    else
                        spikewave = nanmean(cwaves);
                        err_val   = nanstd(cwaves);
                        h_val     = spikewave + err_val;
                        l_val     = spikewave - err_val;
                    end
                    
                    xvec  = 1:length(spikewave);
 
%                     hline(double(plx.SpikeChannels(curr_chan).Threshold)/1000, 'Color', 'black', ...
%                           'LineWidth', 1.5, 'LineStyle',':');
                    plot([xvec(1),xvec(end)], [double(plx.SpikeChannels(curr_chan).Threshold)/1000, double(plx.SpikeChannels(curr_chan).Threshold)/1000]  , ...
                         'Color', 'black', 'LineWidth', 1.5, 'LineStyle',':');
                      
                    patch([xvec, fliplr(xvec)],[h_val fliplr(l_val)],[0.8,0.7,0.7], ...
                           'FaceColor', ptch_col(cunit+1,:), 'EdgeColor', ptch_col(cunit+1,:));
                    phndl(u) = plot(xvec, spikewave, 'k', 'LineWidth', 2, 'Color', main_col(cunit+1,:));
                    plot(xlim,[0,0]);
                    
                    if(minY > min(l_val)) 
                       minY = min(l_val);
                    end
                    if(maxY < max(h_val)) 
                       maxY = max(h_val);
                    end 
                 end
            end
        end
    end
end

% adjust Y-axis to be the same for all channels
for(p=1:length(mplots))
    subplot(mplots(p));
    hold on; box on; axis tight;
    set(gca,'TickDir','out','XTickLabel',[]);
    
    %ofs = 0.05* (maxY - minY);
    %ylim([minY-ofs, maxY+ofs]);
    xlim([xvec(1), xvec(end)]);
    yax = ylim;
    text(0.975*xvec(end), yax(1),  sprintf('DSP%02d%c', chnlst(p)), 'FontSize', 12, 'Interpreter', 'none', ...
        'HorizontalAlignment','right','VerticalAlignment','bottom', 'Rotation',0);
end

% ____________________________________________________________________________ %
%% add information to plots
subplot(mplots(1));
tiHan = title(plxname, 'FontSize', 12,'Interpreter','none');
tiPos = get(tiHan, 'position'); % axis x, y value of title position
xyrange = axis;
set(tiHan, 'position', tiPos + [0 -0.05 * (xyrange(4) - xyrange(3)) 0]); % move title up

% spk_sfx = 'a':'k'; % there should not be much more needed
% 
% if(plot_unsorted == 1) % remove unsorted spike wave forms
%     lgdtext = {unsorted};
%     nunit   = length(phndl)-1;
%     cnt = 1;
% else
%     lgdtext = {};
%     nunit   = length(phndl);
%     cnt = 0;
% end
% for(l=1:nunit)
%     lgdtext(cnt+l) = {spk_sfx(l)};
% end
% legend(phndl,lgdtext, 'Location', 'northeast', 'Box','off', 'Orientation','horizontal')

% ____________________________________________________________________________ %
%% save plot as file
if(~isempty(flnm))
    [~,~,fltyp] = fileparts(flnm);
    fltyp(1)=[];
    
    hgexport(gcf, flnm, hgexport('factorystyle'), 'Format', fltyp);
    
    close(fhndl);
end
