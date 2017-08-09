function [CMABurstInfo,TypeOfSpikeCell] = BurstsAndClusters_BurstPlotterNoRounding(DataCell,ISIcell,CellsInterested);%,BurstAlpha,TailAlpha)

%flags for saving figures and excel files
saveCMAflag = 0 ;   %saving CMA histogram figures as JPEG -set 1 to activate
saveBurstsCMAflag = 0;%saving detection result figures as JPEG -set 1 to activate
xlsCMAwriteflag = 0;%saving burst information on an excel file -set 1 to activate
adaptiveAlpha = 1; %making the CMA detection process adaptive if you do not want to set alpha values yourself. If you want manual alpha then set BurstAlpha and TailAlpha between 0 to 1

% XLSCell = cell(3,CellsInterested);
for i = 1 : length(CellsInterested)
    %% Plotting CMA curve and thresholds
    CMAcurve = ISIcell{CellsInterested(i),9};
    CumulativeHistogram = ISIcell{CellsInterested(i),8};
    ISIHistogram = ISIcell{CellsInterested(i),7};
    if adaptiveAlpha == 1 %for auto alpha selection according to skewness
        Skw = ISIcell{CellsInterested(i),11};
        if Skw < 1
            BurstAlpha=1;
            TailAlpha=0.5;
        elseif Skw < 4
            BurstAlpha=0.7;
            TailAlpha=0.5;
        elseif Skw < 9
            BurstAlpha=0.5;
            TailAlpha=0.3;
        else
            BurstAlpha=0.3;
            TailAlpha=0.1;
        end
    end
    if max(CumulativeHistogram) ~= 0    %ERROR CASE!!!!!!
        MaxCMA =  max(CMAcurve);
        
        MaxPointCMA = max(find(CMAcurve==max(CMAcurve)));
        
        BurstDistants = abs(CMAcurve(MaxPointCMA:end) - (BurstAlpha * max(MaxCMA)));
        BurstThreshold = MaxPointCMA + max(find(BurstDistants==min(BurstDistants))) -1;
        
        if exist('TailAlpha','var') == 1
            TailDistants = abs(CMAcurve(BurstThreshold:end) - (TailAlpha * max(MaxCMA)));
            TailThreshold = BurstThreshold + max(find(TailDistants==min(TailDistants))) -1;
        end
        
        %% Plotting ISI Histogram, Cumulative Histogram, CMA curve and Burts Threshold
        figure(1)
        semilogy(CumulativeHistogram,'b','LineWidth',2)
        hold
        bar(ISIHistogram, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5])
        semilogy(CMAcurve,'r','LineWidth',2)
        %     pause
        line(1:max(CumulativeHistogram)) = BurstThreshold;
        plot(line,1:max(CumulativeHistogram),'k--','LineWidth',2)
        
        line(1:max(CumulativeHistogram)) = TailThreshold;
        plot(line,1:max(CumulativeHistogram),'r--','LineWidth',2)
        axis([0 1000 0 max(CumulativeHistogram)])
        xlabel('ISI(ms)','FontSize',12)
        ylabel('Spike Count','FontSize',12)
        FigureTitle = DataCell{CellsInterested(i),1};
        title(FigureTitle)
        %                         pause
        filename = FigureTitle;
        
        %   Checking the flags for saves
        
        if saveCMAflag == 1
            saveas(figure(1),[filename '_CMA.jpg'])
            saveas(figure(1),[filename '_CMA.fig'])
            %             pause
        end
        
        %         close
        clear line
        %% Plotting Spikes and Bursts
        
        TimeDuration = 700; %Default duration of a recording is assumed as 700 seconds in case there is no information.
        TimeDuration_ms = TimeDuration*1000;
        SpikesFigure = zeros(1,TimeDuration_ms);
        
        temp = round(DataCell{CellsInterested(i),3});
        
        for j = 1: length(temp)
            SpikesFigure(temp(j)) = 1;
        end
        
        figure(2)
        AxisInSeconds = linspace(1,TimeDuration,TimeDuration_ms);
        plot(AxisInSeconds,SpikesFigure)
        axis([1 TimeDuration 0 3 ])
        xlabel('Time(s)','FontSize',12)
        set(gca,'ytick',[])
        title(FigureTitle)
        %         skwn=num2str(ISIcell{i, 11})  added script to plot skwness values
        %                                       on detected burst figures
        %         text(20,3,skwn)
        
        hold
        [TypeOfSpike,BurstStarts,BurstEnds,BurstDurations,NumSpikesInBursts,numberBursts, spikes] = ...
            BurstDetectNoRounding(DataCell, CellsInterested(i), BurstThreshold, TailThreshold);
        
        
        if sum(TypeOfSpike) ~= 0
            
            SpikesIncludedBursts = spikes(find(TypeOfSpike == 1)) ; %Index of the all spikes included to bursts
            SpikesExcludedBursts = spikes(find(TypeOfSpike == 0)) ; %Index of the all spikes excluded from bursts
            BurstLines = zeros(1,length(AxisInSeconds));
            
            if length(SpikesIncludedBursts) ~= 0
                
                for z = 1 : length(BurstStarts)
                    BurstLines(BurstStarts(z):BurstEnds(z)) = 2.3;
                end
                
                plot(AxisInSeconds,((BurstLines./BurstLines).*1.5),'Color','black','LineWidth',2);
                %                                 pause
                %   Checking the flags for saves
                
                if saveBurstsCMAflag == 1
                    saveas(figure(2),[filename '_burstsAlphaAdaptive.jpg'])
                    saveas(figure(2),[filename '_burstsAlphaAdaptive.fig'])
                    %                                 pause
                end
                
                CMABurstInfo{i,1}= FigureTitle;
                CMABurstInfo{i,2} = numberBursts;
                CMABurstInfo{i,3} = sum(BurstDurations) / numberBursts /1000; %in seconds
                CMABurstInfo{i,4} = sum(NumSpikesInBursts)/ numberBursts;
                
            else
                CMABurstInfo{i,1}= FigureTitle;
                CMABurstInfo{i,2} = '-';
                CMABurstInfo{i,3} = '-';
                CMABurstInfo{i,4} = '-';
                %         pause
            end
        else
            CMABurstInfo{i,1}= FigureTitle;
            CMABurstInfo{i,2} = '-';
            CMABurstInfo{i,3} = '-';
            CMABurstInfo{i,4} = '-';
        end
        %         close
        
        
        
        %     index=i;
        %%Sending info to Excel
        %     ExcellBurstInfoRecorder(CellsInterested, bursts,filename,BurstThreshold,index)
        %
        %     %Recording
        %     XLSCell {i,1} = filename;
        %
        %     %Number of bursts
        %     burstspikes = find(bursts==1);
        %     XLSCell {i,2} = length(find(diff(burstspikes)> BurstThreshold));
        %
        %     %bursts per minute
        %     XLSCell {i,3} = XLSCell {i,2}/300;
    else % NO BURST TO DETECT SHOW ONLY SPIKES !!!
        TimeDuration = 700; %Default duration of a recording is assumed as 300 seconds in case there is no information.
        TimeDuration_ms = TimeDuration*1000;
        SpikesFigure = zeros(1,TimeDuration_ms);
        
        temp = round(DataCell{CellsInterested(i),3});
        
        for j = 1: length(temp)
            SpikesFigure(temp(j)) = 1;
        end
        
        figure(2)
        AxisInSeconds = linspace(1,TimeDuration,TimeDuration_ms);
        plot(AxisInSeconds,SpikesFigure)
        axis([1 TimeDuration 0 3 ])
        xlabel('Time(s)','FontSize',12)
        set(gca,'ytick',[])
        FigureTitle = DataCell{CellsInterested(i),1};
        title(FigureTitle)
        %             pause
        %         close
        CMABurstInfo{i,1}= FigureTitle;
        CMABurstInfo{i,2} = '-';
        CMABurstInfo{i,3} = '-';
        CMABurstInfo{i,4} = '-';
        
        spiketimes = DataCell{CellsInterested(i),3};
        TypeOfSpike = zeros(1,length(spiketimes));
        
        if xlsCMAwriteflag == 1
            xlswrite('CMABurstInfoAlphaAdaptive.xls', CMABurstInfo,'Info');
        end
    end
    
    TypeOfSpikeCell{i}=TypeOfSpike;
    
end


