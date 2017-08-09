binwidth = 20;  % for all RT distributions
figure(1); clf;

set(gcf,'Units','Normalized','Position',[0.01 0.03 .98 .89]);

%==================================
%  NSS RT distribution
%==================================

%subplot(3,2,1);
subplot('Position',[0.05 0.67 0.42 0.23]);

edges = []; 
edges = [min(dataSim)-binwidth:binwidth:max(dataSim)+binwidth];
bar(edges, histc(dataSim, edges));
title('NSS RT distribution','FontWeight','Bold')
xlabel('RT (ms)');
ylabel('count');
%%ADDED MLM
%set(gca,'Units','Normalized','Position',[0

%==================================
% plot inhib func
%==================================

%subplot(3,2,2);
subplot('Position',[0.55 0.67 0.42 0.23]);

plot(SSDall, inhibSim, 'k.-', 'linewidth', 3);
hold on;
title('Inhibition function','FontWeight','Bold');
xlabel('SSD (ms)');
ylabel('P(signal-respond)');
ylim([0 1]);
text(80, .9, ['SSRT =' num2str(ceil(overallMeanSSRT))])

%==================================
%  RT distributions - non-canceled
%    and plot CDF plots  
%==================================

%subplot(3,2,3);
subplot('Position',[0.05 0.36 0.42 0.23]);

for x = 1:1:length(SSDall)
    RTdata2 = [];
    RTdata2 = SS_rtSim(x).rt; % NC trials
    
    if ~isempty(RTdata2)
        plotme = [];edges = [];
        [plotme edges] = cdfFunc(RTdata2, 1, min(RTdata2)-10, max(RTdata2));
            switch x
                case 1, plot(edges, plotme, 'r', 'linewidth', 3); hold on;                    
                case 2, plot(edges, plotme, 'y', 'linewidth', 3); hold on;
                case 3, plot(edges, plotme, 'g', 'linewidth', 3); hold on;
                case 4, plot(edges, plotme, 'c', 'linewidth', 3); hold on;
                case 5, plot(edges, plotme, 'b', 'linewidth', 3); hold on;
                case 6, plot(edges, plotme, 'r:', 'linewidth', 3); hold on;
                case 7, plot(edges, plotme, 'g:', 'linewidth', 3); hold on;
            end
    end               
    title('Signal-Respond RT distributions','FontWeight','Bold')
end
NumSSDs = x;

% plot NSS distribution
plotData = [];
[plotData CDFedges] = cdfFunc(dataSim, 1, min(dataSim)-10, max(dataSim));
plot(CDFedges, plotData, 'k', 'linewidth', 3);
ylabel('p(RT< x)')
xlabel('Reaction time (ms)')

switch NumSSDs
    case 1, legend('SSD1', 'NSS', 'location', 'EastOutside');
    case 2, legend('SSD1', 'SSD2', 'NSS', 'location', 'EastOutside');;
    case 3, legend('SSD1', 'SSD2', 'SSD3', 'NSS', 'location', 'EastOutside');
    case 4, legend('SSD1', 'SSD2', 'SSD3', 'SSD4', 'NSS', 'location', 'EastOutside');
    case 5, legend('SSD1', 'SSD2', 'SSD3', 'SSD4', 'SSD5', 'NSS', 'location', 'EastOutside');
    case 6, legend('SSD1', 'SSD2', 'SSD3', 'SSD4', 'SSD5', 'SSD6', 'NSS', 'location', 'EastOutside');
    case 7, legend('SSD1', 'SSD2', 'SSD3', 'SSD4', 'SSD5', 'SSD6', 'SSD7', 'NSS', 'location', 'EastOutside');
end



%==================================
%  cancel time  
%==================================
%subplot(3,2,4)
subplot('Position',[0.55 0.36 0.42 0.23]);
bar([SSDall; SSDall]', [meanGOunitCancelTime; meanSTOPunitCancelTime]')
title('GO and STOP unit cancel times','FontWeight','Bold')
legend('GO unit', 'STOP unit', 'location', 'EastOutside')
ylabel('cancel time')
xlabel('SSD (ms)')

%==================================
%  activation functions - signal respond  
%==================================
ssd = 3;  % plot which SSD?

%subplot(3,2,5); 
subplot('Position',[0.05 0.05 0.42 0.23]);

ind = find(ActFunc(ssd).NC>1000);   % only plot data up to the time that threshold is reached
if isempty(ind) ind=600; end
plot(ActFunc(ssd).NC(1:ind(1)), 'k--', 'linewidth', 3); hold on;
ind = find(ActFunc(ssd).NSSnc>1000);
if isempty(ind) ind=600; end
plot(ActFunc(ssd).NSSnc(1:ind(1)), 'k', 'linewidth', 1); hold on;
plot(ActFunc(ssd).STOPnc, 'k--', 'linewidth', 1); hold on;
plot([SSDall(ssd)+SSRT_by_SSD(ssd) SSDall(ssd)+SSRT_by_SSD(ssd)],[0 1050], 'k')  % plot SSRT
plot([0 length(ActFunc(ssd).Can)],[1000 1000], 'k')  % plot threshold
legend('GOunit_s_i_g_R_e_s_p_o_n_d', 'GOunit_l_a_t_M_a_t_c_h_N_S_S', 'STOPunit', 'SSRT', 'location', 'EastOutside')
ylim([0 1050])
xlim([0 400])
title('Signal Respond Comparison','FontWeight','Bold')
ylabel('GO and STOP activation functions')
xlabel('time relative to trial start (ms)')

%==================================
%  activation functions - signal inhibit  
%==================================

%subplot(3,2,6);
subplot('Position',[0.55 0.05 0.42 0.23]);

ind = find(ActFunc(ssd).Can>1000);   % only plot data up to the time that threshold is reached
if isempty(ind) ind=600; end
plot(ActFunc(ssd).Can(1:ind(1)), 'k', 'linewidth', 3); hold on;
ind = find(ActFunc(ssd).NSSc>1000);   
if isempty(ind) ind=600; end
plot(ActFunc(ssd).NSSc(1:ind(1)), 'k', 'linewidth', 1); hold on;
plot(ActFunc(ssd).STOPc, 'k--', 'linewidth', 1); hold on;
plot([SSDall(ssd)+SSRT_by_SSD(ssd) SSDall(ssd)+SSRT_by_SSD(ssd)],[0 1050], 'k')  % plot SSRT
plot([meanGOunitCancelTime(ssd)+SSDall(ssd)+SSRT_by_SSD(ssd) meanGOunitCancelTime(ssd)+SSDall(ssd)+SSRT_by_SSD(ssd)],...
    [0 1050], 'k:')  % plot cancel time
plot([0 length(ActFunc(ssd).Can)],[1000 1000], 'k')  % plot threshold
legend('GOunit_s_i_g_I_n_h_i_b', 'GOunit_l_a_t_M_a_t_c_h_N_S_S', 'STOPunit', 'SSRT', 'cancel time', 'location', 'EastOutside')
ylim([0 1050])
xlim([0 400])
title('Signal Inhibit Comparison','FontWeight','Bold')
ylabel('GO and STOP activation functions')
xlabel('time relative to trial start (ms)')

H_MMenuPlot = 1;
titlename = 'From: Boucher, L., Palermi, T.J., Logan, G., & Schall, J. D. (2007). Psychological Review.';
% pos_figuretitle=[0 .95 1 .05];
pos_figuretitle=[0 .97 1 .03];
H_Figure_Title =axes('parent',H_MMenuPlot,...
    'position',pos_figuretitle,...
    'box','off',...
    'linewidth', 0.01,...
    'XTick',[],...
    'YTick',[],...
    'color',[1 1 1],...
    'Tag','Figure_Title');
UserFont ='Helvetica';

text(0.5,0.5, titlename,...
    'Fontname',UserFont,...
    'fontsize',[10],...
    'horizontalalignment','center');


% titlename = 'this is the title';
% FigureTitle(titlename)
