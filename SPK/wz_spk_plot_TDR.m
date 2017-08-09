function wz_spk_plot_TDR(TDT, plotwin, rast, fignm)
% ToDo: plot markers for clip times + histograms

%% check arguments and set default
if(exist('plotwin','var') == 0 || isempty(plotwin) == 1)
    plotwin(1) = max([TDT.spk1.timewindow(1); TDT.spk2.timewindow(1)]);
    plotwin(2) = min([TDT.spk1.timewindow(2); TDT.spk2.timewindow(2)]);
end

if(exist('rast','var') == 0 || isempty(rast) == 1)
    rast = 0;
end

% minimum requirement of consecutive significant time bins
if(~exist('fignm','var') || isempty(fignm))
    do_save = 0;
    fignm = [];
else
    do_save = 1;
end

%% set figure properties
% fig = figure('units','normalized','outerposition',[0 0 1 1], 'Renderer', 'Painters', 'Color', [1 1 1]);
fig = figure('Position', [0 0 1600 1000], 'Renderer', 'Painters', 'Color', [1 1 1]);

%% define plot colors
col_sp1_main  = [  9  43   0]./255;
col_sp1_area  = [134 236 107]./255;
col_sp1_2nd   = [ 11  55   0]./255;

col_sp2_main  = [  1   5  45]./255;
col_sp2_area  = [109 121 235]./255;
col_sp2_2nd   = [  1   5  45]./255;

col_diff_main = [ 10  10  10]./255;
col_diff_area = [128 128 128]./255;
col_diff_2nd  = [ 25  25  25]./255;

col_sig_main  = [255 125   0]./255;
col_sig_area  = [255 190  56]./255;
col_sig_2nd   = [255 125   0]./255;

conf_col_area = [255 255 145]./255;
conf_col_line = [255 255   0]./255;

TDTest = TDT.TDT(1);

if(length(TDT.TDT) > 4)
    bwdth = 5; 
    [f,xi] = ksdensity(TDT.TDT, 'width', bwdth); 
    plot(xi,f.*(length(TDT.TDT)*bwdth),'red', 'linewidth',2.5); 

    TDTconfint = prctile(TDT.TDT,[5 95]);
end

%%  raster
if(rast == 1)
    subplot('Position',[0.0375 0.8125 0.6 0.175]);
    hold on;
    maxTrial = TDT.spk1.nTrials+TDT.spk2.nTrials+1;
    ctrial   = maxTrial;

    if(isfield(TDT.spk1,'spiketimes') && isfield(TDT.spk2,'spiketimes') )
        t1 = ctrial;
        for(s=1:TDT.spk1.nTrials)
            ctrial  = ctrial - 1;
            tpos = TDT.spk1.trial_order(s);
            cspikes = TDT.spk1.spiketimes(tpos,:);
            cspikes(isnan(cspikes)) = [];
            if(~isempty(cspikes))
                if(length(TDT.spk1.clip) == TDT.spk1.nTrials)
                    psp = cspikes<=TDT.spk1.clip(tpos);
                    if(any(psp))
                        plot(cspikes(psp),ctrial,'.', 'color', col_sp1_2nd);
                    end
                    psp = cspikes>TDT.spk1.clip(tpos);
                    if(any(psp))
                        plot(cspikes(psp),ctrial,'.', 'color', [0.6 0.6 0.6]);
                    end
                    plot([TDT.spk1.clip(tpos), TDT.spk1.clip(tpos)],[ctrial-0.5, ctrial+0.5],'-r', 'LineWidth',2.5);
                else
                    plot(cspikes,ctrial,'.', 'color', col_sp1_2nd);
                end
            end
        end
        plot([median(TDT.spk1.clip), median(TDT.spk1.clip)],[t1,ctrial],'-r', 'LineWidth',2.5);

        ctrial = ctrial - 2;
        sepln = refline(0,ctrial);
        set(sepln,'Color','red', 'LineWidth',2);
        ctrial = ctrial - 1;

        t2 = ctrial;
        for(s=1:TDT.spk2.nTrials)
            ctrial  = ctrial - 1;
            tpos = TDT.spk2.trial_order(s);
            cspikes = TDT.spk2.spiketimes(tpos,:);
            cspikes(isnan(cspikes)) = [];
            if(~isempty(cspikes))
                if(length(TDT.spk2.clip) == TDT.spk2.nTrials)
                    psp = cspikes<=TDT.spk2.clip(tpos);
                    if(any(psp))
                        plot(cspikes(psp),ctrial,'.', 'color', col_sp2_2nd);
                    end
                    psp = cspikes>TDT.spk2.clip(tpos);
                    if(any(psp))
                        plot(cspikes(psp),ctrial,'.', 'color', [0.6 0.6 0.6]);
                    end
                    plot([TDT.spk2.clip(tpos), TDT.spk2.clip(tpos)],[ctrial-0.5, ctrial+0.5],'-r', 'LineWidth',2.5);
                else
                    plot(cspikes,ctrial,'.', 'color', col_sp2_2nd);
                end
            end
        end
        plot([median(TDT.spk2.clip), median(TDT.spk2.clip)],[t2,ctrial],'-r', 'LineWidth',2.5);
    end

    xlim(plotwin);
    ylim([0, maxTrial]);
%    set(gca,'fontsize',14);
    axis off;
end

%% main plot response estimates
subplot('Position',[0.0375 0.505 0.6 0.305]);
hold on;
    
     m_vals = TDT.spk1.meandensity;
     l_vals = TDT.spk1.meandensity - TDT.spk1.stedensity;
     h_vals = TDT.spk1.meandensity + TDT.spk1.stedensity;
     x_vals = TDT.spk1.xtime;
       
     m_vals2 = TDT.spk2.meandensity;
     l_vals2 = TDT.spk2.meandensity - TDT.spk2.stedensity;
     h_vals2 = TDT.spk2.meandensity + TDT.spk2.stedensity;
     x_vals2 = TDT.spk2.xtime;
     
%     maxpos = min([find(isnan(m_vals),  1, 'first'), find(isnan(l_vals),  1, 'first'),find(isnan(h_vals),  1, 'first'), ...
%                   find(isnan(m_vals2), 1, 'first'), find(isnan(l_vals2), 1, 'first'),find(isnan(h_vals2), 1, 'first')]);
    maxpos = min([find(TDT.spk1.Ndense < 4,  1, 'first'), find(TDT.spk2.Ndense < 4, 1, 'first')]);

    m_vals(maxpos:end) = [];
    l_vals(maxpos:end) = [];
    h_vals(maxpos:end) = [];
    x_vals(maxpos:end) = []; 

    m_vals2(maxpos:end) = [];
    l_vals2(maxpos:end) = [];
    h_vals2(maxpos:end) = [];
    x_vals2(maxpos:end) = []; 

    p1 = x_vals  >= plotwin(1) & x_vals  <= plotwin(2);
    p2 = x_vals2 >= plotwin(1) & x_vals2 <= plotwin(2);
    
    if(length(TDT.TDT)>1)
        maxY = max([h_vals2(p2), h_vals(p1)]);
        patch([TDTconfint, fliplr(TDTconfint)],[0 0 1000 1000], conf_col_area, 'EdgeColor', conf_col_line);
        patch(xi,0.2*maxY*f./max(f),'black');

        t0 = vline(TDT.times([1,end]),'color','g', 'linestyle', '--', 'linewidth',2.5);
        t1 = vline(TDTest,'color','r', 'linewidth',1);
        t2 = vline(TDT.medianTDT,'color','b', 'linewidth',1);
        t3 = vline(TDT.meanTDT,  'color','g', 'linewidth',1);
        t4 = vline(TDT.pTDT(2),  'color','m', 'linewidth',1);
        
        legend([t0 t1 t2 t3 t4],'analysis window', 'raw TDT', 'median boot TDT', 'mean boot TDT', 'p based TDT', 'Location', 'NorthWest');
    else
        maxY = max([h_vals(p1), h_vals2(p2)]);
        t0 = vline(TDT.timeWin,'color','g', 'linestyle', '--', 'linewidth',2.5);
        t1 = vline(TDTest,'color','r', 'linewidth',1);
        legend([t0 t1], 'analysis window', 'TDT', 'Location', 'NorthWest');
    end 

    % second response set
    patch([x_vals2, fliplr(x_vals2)],[l_vals2, fliplr(h_vals2)], col_sp2_area/2.5, 'EdgeColor', col_sp2_2nd/2.5);
    plot(x_vals2, m_vals2, 'Color', col_sp2_main./2.5, 'linewidth', 2.5);
       
    % first response set
    patch([x_vals, fliplr(x_vals)],[l_vals, fliplr(h_vals)], col_sp1_area/2.5, 'EdgeColor', col_sp1_2nd/2.5);
    plot(x_vals, m_vals, 'Color', col_sp1_main/2.5, 'linewidth', 2.5);
        
    m_vals = TDT.spkmean1;
    l_vals = TDT.spkmean1 - TDT.spkste1;
    h_vals = TDT.spkmean1 + TDT.spkste1;
    x_vals = TDT.times;
      
    m_vals2 = TDT.spkmean2;
    l_vals2 = TDT.spkmean2 - TDT.spkste2;
    h_vals2 = TDT.spkmean2 + TDT.spkste2;
    x_vals2 = TDT.times;  
    
     % second response set
    patch([x_vals2, fliplr(x_vals2)],[l_vals2, fliplr(h_vals2)], col_sp2_area, 'EdgeColor', col_sp2_2nd);
    plot(x_vals2, m_vals2, 'Color', col_sp2_main, 'linewidth', 2.5);
       
    % first response set
    patch([x_vals, fliplr(x_vals)],[l_vals, fliplr(h_vals)], col_sp1_area, 'EdgeColor', col_sp1_2nd);
    plot(x_vals, m_vals, 'Color', col_sp1_main, 'linewidth', 2.5);  
    
    ylabel('firing rate [spikes/s]', 'FontSize', 14);
%    xlabel('time [ms]');
    set(gca,'TickDir','out');
    set(gca,'fontsize',12);
    xlim(plotwin);
    ylim([0, maxY])
    set(gca,'layer','top');    

%% plot response difference
subplot('Position',[0.038 0.28 0.6 0.185]);
hold on;
    x_vals = TDT.times;
    m_vals = TDT.meandiff;
    h_vals = TDT.meandiff + TDT.stddiff;
    l_vals = TDT.meandiff - TDT.stddiff;
    maxpos = min([find(isnan(m_vals),  1, 'first'), find(isnan(l_vals),  1, 'first'),find(isnan(h_vals),  1, 'first')]);

    m_vals(maxpos:end) = [];
    l_vals(maxpos:end) = [];
    h_vals(maxpos:end) = [];
    x_vals(maxpos:end) = []; 
   
    if(length(TDT.TDT)>1)
        minval = min(l_vals);
        maxval = max(h_vals);
        patch([TDTconfint, fliplr(TDTconfint)],[minval minval maxval maxval], conf_col_area, 'EdgeColor', conf_col_line);
        patch(xi,minval+(0.2*(maxval-minval)*f./max(f)),'black');
        vline(TDT.medianTDT,'color','b', 'linewidth',1);
        vline(TDT.meanTDT,  'color','g', 'linewidth',1);
        vline(TDT.pTDT(2),  'color','m', 'linewidth',1);

        patch([x_vals, fliplr(x_vals)],[l_vals, fliplr(h_vals)], col_diff_area, 'EdgeColor', col_diff_2nd);
    else
        minval = min(m_vals);
        maxval = max(m_vals);        
    end   
    
    vline(TDTest,'color','r', 'linewidth',1);
    vline(TDT.timeWin,'color','g', 'linestyle', '--', 'linewidth',2.5);
    plot(x_vals, m_vals, 'Color', col_diff_main, 'linewidth', 2.5);

    hline(0,'Linewidth',1)
    
%    xlabel('time [ms]');
    ylabel('response difference [spikes/s]', 'FontSize', 14);
    set(gca,'TickDir','out');
    set(gca,'fontsize',12);
    xlim(plotwin);
    ylim([minval maxval]);
    set(gca, 'XTick', []);
    set(gca,'layer','top'); 
    box on 

%% plot p-value TS
subplot('Position',[0.0375 0.075 0.6 0.185]);
hold on;

    x_vals = TDT.times;
    m_vals = log(median(TDT.Pmat,1));
    h_vals = log(prctile(TDT.Pmat,75,1));
    l_vals = log(prctile(TDT.Pmat,25,1));

%     m_vals(maxpos:end) = [];
%     l_vals(maxpos:end) = [];
%     h_vals(maxpos:end) = [];
%     x_vals(maxpos:end) = []; 

    if(length(TDT.TDT)>1)
        minval = min(l_vals);
        maxval = max(h_vals);
        patch([TDTconfint, fliplr(TDTconfint)],[minval minval maxval maxval], conf_col_area, 'EdgeColor', conf_col_line);
        patch(xi,minval+(0.2*(maxval-minval)*f./max(f)),'black');
        vline(TDT.medianTDT,'color','b', 'linewidth',1);
        vline(TDT.meanTDT,  'color','g', 'linewidth',1);
        vline(TDT.pTDT(2),  'color','m', 'linewidth',1);
        
        patch([x_vals, fliplr(x_vals)],[l_vals, fliplr(h_vals)], col_sig_area, 'EdgeColor', col_sig_2nd);
    else
        minval = min(m_vals);
        maxval = max(m_vals);        
    end
    
    vline(TDTest,'color','r', 'linewidth',1);
    vline(TDT.timeWin,'color','g', 'linestyle', '--', 'linewidth',2.5);
    plot(x_vals, m_vals, 'Color', col_sig_main, 'linewidth', 2.5);
    
    p05  = hline(log(0.05), 'color','red','LineStyle',':');
    p01  = hline(log(0.01), 'color','red','LineStyle','--');
    p001 = hline(log(0.001),'color','red','LineStyle','-');
    legend([p05 p01 p001], 'p = 0.05', 'p = 0.01', 'p = 0.001', 'Location', 'SouthEast');

    ylabel('log(p)', 'FontSize', 14);
    xlabel('time [ms]', 'FontSize', 14);
    set(gca,'TickDir','out');
    set(gca,'fontsize',12);
    xlim(plotwin);
    ylim([minval maxval]);
    set(gca,'layer','top'); 
    box on
    
%% plot histogram
if(length(TDT.TDT) > 4)
    subplot('Position',[0.675 0.505 0.3 0.305]);
    hold on;

        minTDT = min(TDT.TDT)-2.5;
        maxTDT = max(TDT.TDT)+2.5;

        Tbin = minTDT-1.5*bwdth : bwdth : maxTDT+1.5*bwdth;

        patch([TDTconfint, fliplr(TDTconfint)],[0 0 length(TDT.TDT) length(TDT.TDT)],conf_col_area, 'EdgeColor', conf_col_line);

        TDThist = histc(TDT.TDT,Tbin);
        bar(Tbin, TDThist, 'k');
        [f,xi,bw] = ksdensity(TDT.TDT, 'bandwidth', bwdth); 
        f = f.*bw*sum(~isnan(TDT.TDT));
        plot(xi,f,'red', 'linewidth',2.5); 

        vline(TDTest,'color','r', 'linewidth',1);
        if(length(TDT.TDT)>1)
            vline(TDT.medianTDT,'color','b', 'linewidth',1);
            vline(TDT.meanTDT,  'color','g', 'linewidth',1);
            vline(TDT.pTDT(2),  'color','m', 'linewidth',1);
        end

        ylabel('count', 'FontSize', 14);
        xlabel('TDT [ms]', 'FontSize', 14);
        xlim([min(Tbin), max(Tbin)]);
        ylim([0 max([max(TDThist), max(f)])]);
        set(gca,'fontsize',12);
        set(gca,'layer','top'); 
end

%% write data
subplot('Position',[0.675 0.075 0.3 0.38]);
hold on;
axis off
if(~isempty(fignm))
    [~, flnm] = fileparts(fignm);
    title(flnm);
end
    


if(length(TDT.TDT) == 1)
    smry_str = {['TDT = ', int2str(TDT.TDT),' ms']};
else
    smry_str = {['Nboot = ', int2str(length(TDT.TDT))] ...
                ['success rate = ', num2str(100*sum(~isnan(TDT.TDT))/length(TDT.TDT),'%.2f'), '%'] ...
                [''],...
                ['raw TDT = ', int2str(TDT.TDT(1)),' ms'], ...
                [''],...
                ['mode boot TDT = ',   int2str(mode(TDT.TDT)),' ms'], ...
                ['median boot TDT = ', int2str(TDT.medianTDT),' ms'], ...
                ['IQR boot TDT = ',    int2str(iqr(TDT.TDT)), ' ms'], ...
                ['P05 boot TDT = ',    int2str(prctile(TDT.TDT, 5)),' ms'], ...
                ['P25 boot TDT = ',    int2str(prctile(TDT.TDT,25)),' ms'], ...
                ['P75 boot TDT = ',    int2str(prctile(TDT.TDT,75)),' ms'], ...
                ['P95 boot TDT = ',    int2str(prctile(TDT.TDT,95)),' ms'], ...
                ['MAD boot TDT = ',    int2str(mad(TDT.TDT,1)),     ' ms'], ...
                [''],...
                ['median TDT p 0.05 = ',  int2str(TDT.medTDTp05), ' ms'], ...
                ['median TDT p 0.001 = ', int2str(TDT.medTDTp001),' ms'], ...
                ['median p distance = ',  int2str(TDT.medpdist),  ' ms'], ...
                [''],...
                ['mean boot TDT = ', num2str(TDT.meanTDT,'%.2f'),    ' ms'], ...
                ['std boot TDT = ',  num2str(nanstd(TDT.TDT),'%.2f'),' ms'], ...
                [''],...
                ['median p based TDT = ', int2str(TDT.pTDT(2)),' ms'] ...
                ['lower p based TDT = ',  int2str(TDT.pTDT(1)),' ms'] ...
                ['upper p based TDT = ',  int2str(TDT.pTDT(3)),' ms'] ...
                };
end

text(0, 1, 1, smry_str, 'FontSize', 14, 'VerticalAlignment', 'top');

%% save figure
if(do_save == 1)
%    saveas(fig, fignm, 'eps');
%     saveas(fig, fignm, 'png');
    
    hgexport(fig, fignm,  ...
    hgexport('factorystyle'), 'Format', 'png'); 
%     hgexport('factorystyle'), 'Format', 'eps'); 
%     hgexport(fig, fignm,  ...
    close all;
end
