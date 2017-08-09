
% *********************************************************************************************
% *********************************************************************************************
% *******   This code was authored by Michael Hill (buffalohill@me.com), it may not     *******   
% *******   be used as a whole or in parts without the written permission of Michael    *******  
% *******   Hill and any use of this code (in whole or in parts) must reference the     ******* 
% *******   paper in which the h-coefficient method was first introduced (contact me    ******* 
% *******   to inquire about appropriate referencing formats for this code).             ****** 
% *********************************************************************************************
% *********************************************************************************************


% The following parameters need to be provided to model and analyze a
% neurone. You can change these but not all combinations of all types of
% parameters may work. To be on the safe side only change one parameter at a
% time and see what effect it has and how stable the script still is.
% Within these limitations you should still be able to model almost any
% kind of neuronal response within reason. To get started try running the 
% code multiple times adding the following changes sequentially and see how
% the results change accordingly (at some point also change doplot to 0 to 
% speed things up):
%
% (1) try changing: noisefreq   to   [0,0];     ->      No noise (homogenous Poisson process)
% (2) try changing: amp         to   3;         ->      High response amplitude
% (2) try changing: mfr         to   30;        ->      High baseline firing rate
% (3) try changing: percent     to   0.66;      ->      Only 66% of all trials evoked a response in this neurone
% (4) try changing: percent     to   -0.66;     ->      Only 66% of all trials evoked a response in this neurone, the remaining 33% evoked inhibition

% Enter here the path to where you want to save the results figures
saveto_path = '//Users/Portahill/Desktop';

%% (1) Set the parameters for you modelled data

plotrange   =   [-1500 2000];   % The range around your TTL (t = 0), which you want to plot in your PSTH in milliseconds
binwidth    =   200;            % The bindwidth of the comparison histogram in your plot in milliseconds
overlap     =   1;              % Plot the histogram with 50% overlapping bins (1) or with non-overlapping bins (0) 
duration    =   30;             % The duration of your whole experiment in minutes (i.e. the length of the raw data)
nstims      =   30;             % Then umber of stimuli (i.e. the number of 'lines' in your raster plots)
percent     =   1;              % The percentage of trials during which the modelled neurone responded to the stimuli (e.g. 1 = 100%, i.e. the neurone responds every time, 0.5 = 50%, i.e. the neurone only fired in response to half of the presented stimuli). You can also enter negative values here, in which case the neurone is inhibited during those trials in which it did not responde (inhibited to the same amount as it's firing rate is increase during responding trials with an obvious floor at 0Hz) 
noisefreq   =   [1,0.1];        % The frequency range of a sine wave, which will be convoluted with your data as a noise signal. The first value is the frequency of the noise, the second value is the amount by which the frequency is randomly  modulated in each trial(i.e. +/- that amount).
noiseamp    =   [1,0.1];        % The amplitude of the above noise sine wave, the first value is the  amplitude of the noise, the second value is the amount by which it is randomly modulated in each trial (i.e. +/- that amount, the actual values of the noiseamp are entered in multiples of the baseline value, see below).
xdata       =   (1:1:10000);    % The sampling frequence (please do not change as the script is not very flexible in this respect)
sigma       =   300;            % The width of the neurons response to a stimulus in ms (i.e. the standard deviation of your gaussian). sigma is entered as a scalar in units of milliseconds, e.g. sigma = 300 (ms). 
rspdealy    =   5450;           % The location of the neuron's response peak in ms relative to the beginning of your xdata range (t = 0 is actually at 5000 here, i.e. in the middle of you xdata range).
amp         =   1.25;           % The amplitude of your response relative to the neurons baseline firing rate. For example: amp = 0   The neuron responds with total inhibition, i.e. 0 Hz amp = 1   The neuron shows no response, i.e. it keeps on firing at baseline amp = 2 	The neuron shows a response with an amplitude double of your baseline firing rate amp is entered as a scalar normalized to baseline, e.g. amp = 2.
bslrange    =   [0,4800];       % The beginning and end of your baseline range (see help expsim for details).
mfr         =   3;              % mfr is the mean firing rate during your baseline range in Hz.
doplot      =   1;              % Set to 1 to see the Monte Carlo simulation in action, set to 0 to increase speed significanly
fast        =   1;              % Mor modelled data such as is the case here, the fast h-coefficient can always be used, for real data this needs to be decided based on the available computational power and the researchers own preferece.
rsp         =   [200,1000];     % The range within which a response is expected.     

%% (2) Generate some Monte Carlo simulated neuronal responses according to the above parameters

[peristim,peristim_cntr,grndtruth,grndtruth_cntr,rawdata] = expsim(duration,nstims,percent,noisefreq,noiseamp,xdata,sigma,rspdealy,amp,bslrange,mfr,doplot);

%% (3) Now plot your results (one PSTH with a response in it and one as a control)

ylims = [-0.1,1.5*amp+1];
semcol = 0.8;
for cnd = 1:2
    switch cnd
        case 1
            grnd = grndtruth;
            peris = peristim;
        case 2
            grnd = grndtruth_cntr;
            peris = peristim_cntr;
    end
    close all
    figure
    % Subplot 1
    subplot(3,1,1);
    hold on
    [kern,hcoeffs,hcoeffs2D] = hcoeff(peris,rawdata,mfr,nstims,rsp,fast);
    if percent == 1
        for a = 1:size(grnd,1)
            plot((xdata-5000)/1000,grnd(a,:),'-','Color',[0.6 0.6 0.6],'LineWidth',1);
        end
    else
        for a = 1:size(grnd,1)-(ceil(abs(percent) * nstims))
            plot((xdata-5000)/1000,grnd(a,:),'-g','LineWidth',1);
        end
        for a = (floor(abs(percent) * nstims))+1:size(grnd,1)
            plot((xdata-5000)/1000,grnd(a,:),'-r','LineWidth',1);
        end
    end
    jbfill((xdata-5000)/1000,nanmean(grnd,1)+(nanstd(grnd,1)/sqrt(size(grnd,1))),nanmean(grnd,1)-(nanstd(grnd,1)/sqrt(size(grnd,1))),[semcol semcol semcol],[semcol-0.1 semcol-0.1 semcol-0.1],[],0.75);
    plot((xdata-5000)/1000,nanmean(grnd,1),'-k','LineWidth',2);
    box off
    axis tight
    grid on
    set(gca,'FontSize',9);
    set(gca,'XColor','w');
    set(gca,'XTickLabel',{});
    xlim([-1.5 2.5])
    ylabel('Ground truth')
    ylim(ylims);
    xlim(plotrange/1000)
    if cnd == 2
        title([{['Control Monte Carlo simulation without a response and a mean firing rate of ',num2str(mfr),' Hz.']}]);
    elseif percent == 1
        title([{['Monte Carlo simulation with ',num2str(abs(percent)*100),'% positive responses with a mean firing rate of ',num2str(mfr),' Hz,']};...
            {[' response center at ',num2str(((rspdealy-5000)/1000)),' sec, sigma of ',num2str(sigma/1000),' sec and an amplitude of ',num2str(abs(amp)+1),' times baseline.']}]);
    elseif percent > 0
        title([{['Monte Carlo simulation with ',num2str(abs(percent)*100),'% positive responses with a mean firing rate of ',num2str(mfr),' Hz,']};...
            {[' response center at ',num2str(((rspdealy-5000)/1000)),' sec, sigma of ',num2str(sigma/1000),' sec and an amplitude of ',num2str(abs(amp)+1),' times baseline.']};...
            {['Green = positive responses, red = remaining (flat) responses.']}]);
    elseif percent < 0
        title([{['Monte Carlo simulation with ',num2str(abs(percent)*100),'% positive responses with a mean firing rate of ',num2str(mfr),' Hz,']};...
            {[' response center at ',num2str(((rspdealy-5000)/1000)),' sec, sigma of ',num2str(sigma/1000),' sec and an amplitude of ',num2str(abs(amp)+1),' times baseline.']};...
            {['Green = positive responses, red = remaining responses with same negative parameters']}]);
    end
    loc1 = subplot(3,1,1);
    set(loc1,'position',[0.1300,0.6254,0.7750,0.2996]);
    % Subplot 2
    subplot(312)
    prelimbins = -(length(xdata)/2):binwidth/2:(length(xdata)/2);     % create a sequence of binwidth of half width for your histogram with overlapping binwidth (let the histogram span your whole (length(xdata)/2))
    h = 0;
    for a = 1:size(peris,1)
        % get a normal histogram of your spikes
        histo(a,:) = hist(peris{a,1},(-((length(xdata)/2)-(binwidth/2)):binwidth:((length(xdata)/2)-(binwidth/2))));
        % at the same time also get a histogram with overlapping binwidth (the binwidth remain the same size as above but overlab by half their width):
        clear prelimhisto prelimovhist
        prelimhisto = histc(peris{a,1},prelimbins);
        prelimovhist = prelimhisto + [prelimhisto(2:end),0];
        movhisto(a,:) = prelimovhist(1,1:end-2);                     % get rid of the last two entries as they are edge effects (see help histc for details)
        for l=1:size(peris{a,1},2)
            if peris{a,1}(1,l) >= plotrange(1,1) && peris{a,1}(1,l) <= plotrange(1,2)
                line([peris{a,1}(1,l) peris{a,1}(1,l)],[h h+0.6],'color','k','LineWidth',2);              % draw a black line
            end
        end
        h = h+1;
    end
    set(gca,'XColor','w');
    set(gca,'XTickLabel',{});
    ylabel('Rasterplot')
    xlim(plotrange)
    % Subplot 3
    subplot(3,1,3);
    hold on
    if overlap == 1  
        xplot = (-((length(xdata)/2)-(binwidth/2)):binwidth/2:((length(xdata)/2)-(binwidth/2)))./1000;
        bar(xplot,(nanmean(movhisto*(1000/binwidth)/mfr)),'FaceColor',[0.6 0.6 0.6],'EdgeColor','none');
    else
        xplot = (-((length(xdata)/2)-(binwidth/2)):binwidth:((length(xdata)/2)-(binwidth/2)))./1000;
        bar(xplot,nanmean(histo*(1000/binwidth)/mfr,1),'FaceColor',[0.6 0.6 0.6],'EdgeColor','none');
    end
    hold on
    plot(kern(1,:),(kern(2,:)/mfr),'-k','LineWidth',2);
     xlabel([{['Measured h-coefficients (a+b)/c:']};...
            {['h(max) = (',num2str(hcoeffs2D(1,2)),'+',num2str(hcoeffs2D(1,1)),')/',num2str(hcoeffs2D(1,3)),' = ',num2str(roundn(hcoeffs(1,1),-2)),',']};...
            {['h(min) = (',num2str(hcoeffs2D(2,2)),'+',num2str(hcoeffs2D(2,1)),')/',num2str(hcoeffs2D(2,3)),' = ',num2str(roundn(hcoeffs(1,2),-2)),'.']}]);
    ylabel('PSTH')
    ylim([0,max(nanmean(histo*(1000/binwidth)/mfr,1))+0.5]);
    xlim(plotrange/1000)
    origloc1 = subplot(312);
    origloc2 = subplot(313);
    loc1 = get(origloc1,'position');
    loc2 = get(origloc2,'position');
    loc3 = [loc1(1,1) loc2(1,2)+loc2(1,4) loc1(1,3) loc1(1,4)+(loc1(1,2)-(loc2(1,2)+loc2(1,4)))];
    set(origloc1,'position',loc3);
    set(gcf,'Position',[32,241,560,838],'PaperPositionMode','auto','color','w')
    cl = clock;
    prf1 = sprintf('%s-%d-%d-%d',date,cl(4),cl(5),round(cl(6)));
    if cnd == 1
        cd(saveto_path);
        saveas(gcf,char(sprintf('RspCntr-%d_Nstims-%d_mfr-%d_sigma-%d_amp-%s_percent-%d_Nfreq-%s-%s_Namp-%s-%s___%s.fig',round(rspdealy),round(nstims),round(mfr),round(sigma),num2str(amp),round(percent * 100),num2str(roundn(noisefreq(1,1),-1)),num2str(roundn(noisefreq(1,2),-1)),num2str(roundn(noiseamp(1,1),-1)),num2str(roundn(noiseamp(1,2),-1)),prf1)),'fig');
        cd('/Users/Portahill/Desktop');
        print(char(sprintf('RspCntr-%d_Nstims-%d_mfr-%d_sigma-%d_amp-%s_percent-%d_Nfreq-%s-%s_Namp-%s-%s___%s.jpg',round(rspdealy),round(nstims),round(mfr),round(sigma),num2str(amp),round(percent * 100),num2str(roundn(noisefreq(1,1),-1)),num2str(roundn(noisefreq(1,2),-1)),num2str(roundn(noiseamp(1,1),-1)),num2str(roundn(noiseamp(1,2),-1)),prf1)),'-djpeg','-r200','-zbuffer')
        close all
    else
        cd(saveto_path);
        saveas(gcf,char(sprintf('ContCntr-%d_Nstims-%d_mfr-%d_sigma-%d_amp-%s_percent-%d_Nfreq-%s-%s_Namp-%s-%s___%s.fig',round(rspdealy),round(nstims),round(mfr),round(sigma),num2str(amp),round(percent * 100),num2str(roundn(noisefreq(1,1),-1)),num2str(roundn(noisefreq(1,2),-1)),num2str(roundn(noiseamp(1,1),-1)),num2str(roundn(noiseamp(1,2),-1)),prf1)),'fig');
        cd('/Users/Portahill/Desktop');
        print(char(sprintf('ContCntr-%d_Nstims-%d_mfr-%d_sigma-%d_amp-%s_percent-%d_Nfreq-%s-%s_Namp-%s-%s___%s.jpg',round(rspdealy),round(nstims),round(mfr),round(sigma),num2str(amp),round(percent * 100),num2str(roundn(noisefreq(1,1),-1)),num2str(roundn(noisefreq(1,2),-1)),num2str(roundn(noiseamp(1,1),-1)),num2str(roundn(noiseamp(1,2),-1)),prf1)),'-djpeg','-r200','-zbuffer')
    end
    close all
end