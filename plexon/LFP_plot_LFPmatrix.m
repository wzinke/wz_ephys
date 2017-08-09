function [fhndl, stim_plots, sacc_plots] = LFP_plot_LFPmatrix(SessDir, odir, oext)
% Plot LFP signals for all channels of an U-probe.
%  
% DESCRIPTION 
%   The LFP data of all 24 channels of an U-probe will be plotted and aligned to
%   stimulus onset or saccade onset. Currently, target location is split into
%   left (red) versus right (blue).
%
%
% SYNTAX 
% 
%   LFP_plot_LFPmatrix(SessDir, odir, oext)
%
%   Input:
%
%       SessDir         path to session directory which has to contain an 'LFP' subdirectory
%
%       odir            Output directory to save the summary plots (default is <SessDir>)
%
%       oext            extension for the file format to save the figure. If not
%                       specified no data will be written
%
%
%   Output:
%
%       fhndl           handle for the figure
%
%       stim_plots      vector with handles to all subplots with stimulus aligned data
%
%       sacc_plots      vector with handles to all subplots with saccade aligned data
%
% ......................................................................... 
% wolf zinke, wolfzinke@gmail.com 
%
% $Created : 15-Feb-2015 by wolf zinke

stim_win = [-150, 650];
sacc_win = [-650, 150];

% if not specified use GUI to get the file
if(~exist('SessDir','var') || isempty(SessDir))
    SessDir = uigetdir('path to session data');
end

if(~exist('odir','var') || isempty(odir))
    odir = SessDir;
end

if(~exist('oext','var'))
    oext = [];
else
    pmkdir(odir);
end

% ____________________________________________________________________________ %
%% get extension list
[~, cSess] = fileparts(SessDir);
cLFP = sprintf('LFP%02d', 1);

lst = dir(fullfile(SessDir,'LFP',cLFP,[cSess,'_',cLFP,'_*.mat']));

num_fls = length(lst);

for(i=1:num_fls)
    [~,cdata] = fileparts(lst(i).name);
    splpos  = find(cdata == '_');
    fileSFX = cdata(splpos(2)+1:end);
    % ____________________________________________________________________________ %
    %% create plot outline:  Stim | Sacc  || Stim | Sacc  ||  Stim | Sacc  ||  Stim | Sacc
    num_rows = 6;
    num_cols = 4;
    figsz    = [0 0 num_cols*400 num_rows*200];
    
    Ystart = linspace(98,3,num_rows+1)./100;
    Ystart(1)   = [];
    
    Xstart = linspace(3,98,num_cols+1)./100;
    Xstart(end) = [];
    
    xwd = min(diff(Xstart))-0.004;
    ywd = min(abs(diff(Ystart)))-0.006;
    
    plwdth=0.475*xwd;
    
    stim_plots = [];
    sacc_plots = [];
    
    fhndl  = figure('Name', [SessDir, '  -  ', fileSFX], 'Position', figsz, 'Renderer', 'Painters');
    
    for(x=1:num_cols)
        cX = Xstart(x);
        for(y=1:num_rows)
            cY = Ystart(y);
            
            % stimulus aligned plot
            stimhdl = subplot('Position',[cX cY plwdth ywd]);
            hold on;
            vline(0,'LineWidth',0.8);
            if(x>1)
                set(gca,'YTickLabel',[]);
            end
            if(y<num_rows)
                set(gca,'XTickLabel',[]);
            end
            
            hold on; box on; axis tight;
            set(gca,'TickDir','out');
            stim_plots = [stim_plots; stimhdl];
            
            % saccade aligned plot
            sacchdl = subplot('Position',[cX+1.025*plwdth cY plwdth ywd]);
            hold on;
            vline(0,'LineWidth',0.8);
            set(gca,'YTickLabel',[]);
            if(y<num_rows)
                set(gca,'XTickLabel',[]);
            end
            
            hold on; box on; axis tight;
            set(gca,'TickDir','out');
            
            sacc_plots = [sacc_plots; sacchdl];
        end
    end
    
    % ____________________________________________________________________________ %
    %% read in LFP channels and plot data
    STimRightLFPmat = [];
    STimLeftLFPmat  = [];
    SaccRightLFPmat = [];
    SaccLeftLFPmat  = [];
    
    maxY = -inf; minY =  inf;
    
    for(c=1:24)
        cLFP  = sprintf('LFP%02d', c);
       
        cfl = fullfile(SessDir,'LFP',cLFP,[cSess,'_',cLFP,'_',fileSFX,'.mat']);
        
        if(~exist(cfl,'file'))
            continue
        end
        
        dt = load(cfl);
        p = find(dt.Task.Correct == 1);

        nTrials = length(p);
        
        lp = dt.Task.TargetLoc(p) == 180;
        rp = dt.Task.TargetLoc(p) ==   0;
        
        if(c==1)
            pltrng = prctile(dt.Task.SRT(p),90);
            medsac = prctile(dt.Task.SRT(p),50);
            
            xstim = stim_win(1):stim_win(2);
            xsacc = sacc_win(1):sacc_win(2);
            
            plstim = xstim > -50       & xstim < pltrng;
            plsacc = xsacc > -1*pltrng & xsacc < 50;
            
            Nr = sum(rp);
            Nl = sum(lp);
        end
                
        xvec = dt.timevec;
        
        HitLFP  = dt.LFP(p,:);
        
        clip_lfp = bsxfun(@ge,xvec,dt.Task.SRT(p));
        
        clipLFP = HitLFP;
        clipLFP(clip_lfp) = NaN;
        
        % allign stimulus and saccade onset
        StimLFP      = nan(nTrials, diff(stim_win)+1);
        StimLFP_clip = nan(nTrials, diff(stim_win)+1);
        SaccLFP      = nan(nTrials, diff(stim_win)+1);
        
        for(t=1:nTrials)
            stim_P = xvec >= stim_win(1) & xvec <= stim_win(2);
            
            StimLFP(t,1:sum(stim_P))      = HitLFP(t,stim_P);
            StimLFP_clip(t,1:sum(stim_P)) = clipLFP(t,stim_P);
            
            sacc_P = xvec-dt.Task.SRT(p(t)) >= sacc_win(1) & xvec-dt.Task.SRT(p(t)) <= sacc_win(2);
            SaccLFP(t,1:sum(sacc_P))      = HitLFP(t,sacc_P);
        end
        
        %% keep LFP in a matrix
        STimRightLFPmat = [STimRightLFPmat; nanmean(StimLFP_clip(rp,plstim))];
        STimLeftLFPmat  = [STimLeftLFPmat; nanmean(StimLFP_clip(lp,plstim))];
        SaccRightLFPmat = [SaccRightLFPmat; nanmean(SaccLFP(rp,plstim))];
        SaccLeftLFPmat  = [SaccLeftLFPmat; nanmean(SaccLFP(lp,plstim))];
        
        
       %% stimulus onset
        subplot(stim_plots(c));
        
        % target left
        [minY, maxY] = plot_lfp(xstim(plstim), StimLFP_clip(lp,plstim), minY, maxY, 'r');
                
        % target right
        [minY, maxY] = plot_lfp(xstim(plstim), StimLFP_clip(rp,plstim), minY, maxY, 'b');
                
        vline(medsac,'color','g');
        xlim([-50; pltrng]);
        
       %% saccade onset
        subplot(sacc_plots(c));        
        % target left
        [minY, maxY] = plot_lfp(xsacc(plsacc), SaccLFP(lp,plsacc), minY, maxY, 'r');
        
        % target right
        [minY, maxY] = plot_lfp(xsacc(plsacc), SaccLFP(rp,plsacc), minY, maxY, 'b');
               
        vline(-1*medsac,'color','g');
        xlim([-1*pltrng, 50]);
    end
    
    for(c=1:24)
        subplot(stim_plots(c));
        ylim([minY, maxY]);
        if(mod(c,num_rows)==1)
            tiHan = title('Stimulus', 'FontSize', 12,'Interpreter','none');
            tiPos = get(tiHan, 'position'); % axis x, y value of title position
            xyrange = axis;
            set(tiHan, 'position', tiPos + [0 -0.05 * (xyrange(4) - xyrange(3)) 0]); % move title up
        end
        Xrng = xlim;
        text(Xrng(2)-0.05*diff(Xrng), minY+(0.01*(maxY-minY)),  sprintf('LFP%02d', c), 'FontSize', 10, 'Interpreter', 'none', ...
            'HorizontalAlignment','right','VerticalAlignment','bottom', 'Rotation',0);

        if(c==1)
            Xrng = xlim;
            text(Xrng(1)+0.05*diff(Xrng), maxY,  [cSess,'  -  ',fileSFX], 'FontSize', 10, 'Interpreter', 'none', ...
                'HorizontalAlignment','left','VerticalAlignment','top', 'Rotation',0);
        end
        
        subplot(sacc_plots(c));
        ylim([minY, maxY]);
        if(mod(c,num_rows)==1)
            tiHan = title('Saccade', 'FontSize', 12,'Interpreter','none');
            tiPos = get(tiHan, 'position'); % axis x, y value of title position
            xyrange = axis;
            set(tiHan, 'position', tiPos + [0 -0.05 * (xyrange(4) - xyrange(3)) 0]); % move title up
        end
        
        if(c==24)
            Xrng = xlim;
            text(Xrng(2)-0.05*diff(Xrng), minY,  ['N(left): ',num2str(Nl),'  -  N(right): ',num2str(Nr)], 'FontSize', 10, 'Interpreter', 'none', ...
                'HorizontalAlignment','right','VerticalAlignment','bottom', 'Rotation',0);
        end
    end
    
    if(~isempty(oext))
        if(strmatch(oext,'svg'))
            plot2svg(fullfile(odir,[cSess,'_',fileSFX,'.',oext]), fhndl);
        else
            hgexport(gcf, fullfile(odir,[cSess,'_',fileSFX,'.',oext]), hgexport('factorystyle'), 'Format', oext);
        end
        close(fhndl);
    end
    
    figsz = [0 0 4*400 1000];
    fhndl  = figure('Name', [SessDir, '  -  ', fileSFX], 'Position', figsz, 'Renderer', 'Painters');

    maxLFP = max(abs([STimRightLFPmat(:);STimLeftLFPmat(:);SaccRightLFPmat(:);SaccLeftLFPmat(:)]));
    
    STimRightLFPmat = STimRightLFPmat ./ maxLFP;
    STimLeftLFPmat  = STimLeftLFPmat  ./ maxLFP;
    SaccRightLFPmat = SaccRightLFPmat ./ maxLFP;
    SaccLeftLFPmat  = SaccLeftLFPmat  ./ maxLFP;
    
    subplot(1,4,1);
    hold on; axis tight;
    set(gca,'TickDir','out');
    
    for(c=1:24)
        plot(xstim(plstim),STimRightLFPmat(c,:)+25-c,'b','LineWidth',2);
        plot(xstim(plstim),STimLeftLFPmat(c,:) +25-c,'r','LineWidth',2);
    end
    vline(0,'color','k');
    vline(medsac,'color','g');
    xlim([-50; pltrng]);
    ylim([0,25])
    title('Stimulus');
    ylabel('Channel');
    xlabel('time [ms]');
    
    subplot(1,4,2);
    hold on; axis tight;
    set(gca,'TickDir','out');
    StimDiffMat = STimLeftLFPmat - STimRightLFPmat;
    StimDiffMat = StimDiffMat ./ max(abs(StimDiffMat(:)));
    
    for(c=1:24)
        plot(xstim(plstim),StimDiffMat(c,:)+25-c,'k','LineWidth',2);
    end
    vline(0,'color','k');
    vline(medsac,'color','g');
    xlim([-50; pltrng]);
    ylim([0,25])
    title('Stimulus: contra vs. ipsi');
    xlabel('time [ms]');
    
    subplot(1,4,3);
    hold on; axis tight;
    set(gca,'TickDir','out');
    
    for(c=1:24)
        plot(xsacc(plsacc),SaccRightLFPmat(c,:)+25-c,'b','LineWidth',2);
        plot(xsacc(plsacc),SaccLeftLFPmat(c,:) +25-c,'r','LineWidth',2);
  
    end
    vline(0,'color','k');
    vline(-1*medsac,'color','g');
    xlim([-1*pltrng, 50]);
    ylim([0,25])
    title('Saccade');
    xlabel('time [ms]');

    subplot(1,4,4);
    hold on; axis tight;
    set(gca,'TickDir','out');
    SaccDiffMat = SaccLeftLFPmat - SaccRightLFPmat;
    SaccDiffMat = SaccDiffMat ./ max(abs(SaccDiffMat(:)));
    
    for(c=1:24)
        plot(xsacc(plsacc),SaccDiffMat(c,:)+25-c,'k','LineWidth',2);
  
    end
    vline(0,'color','k');
    vline(-1*medsac,'color','g');
    xlim([-1*pltrng, 50]);
    ylim([0,25])
    title('Saccade: contra vs. ipsi');
    xlabel('time [ms]');
    
    
    if(~isempty(oext))
        if(strmatch(oext,'svg'))
            plot2svg(fullfile(odir,[cSess,'_',fileSFX,'_arr.',oext]), fhndl);
        else
            hgexport(gcf, fullfile(odir,[cSess,'_',fileSFX,'_arr.',oext]), hgexport('factorystyle'), 'Format', oext);
        end
        close(fhndl);
    end
    
end


function [min_y, max_y] = plot_lfp(xvec, LFP, min_y, max_y, col)

    Rmn  = nanmean(LFP);
    Rste = nanstd(LFP) ./ sqrt(sum(isfinite(LFP)));

    if(max_y < max(Rmn+Rste))
       max_y = max(Rmn+Rste);
    end

    if(min_y > min(Rmn-Rste))
       min_y = min(Rmn-Rste);
    end       

    shadedErrorBar(xvec,Rmn,Rste,col,0)


