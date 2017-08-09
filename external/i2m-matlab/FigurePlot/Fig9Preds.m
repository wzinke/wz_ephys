function Fig9Preds(PexpTab,PthTab,PthNoJ1Tab,PthNoJTab,FeaturesTab,BinIndex,BinSize,TempTab,state,N_r)

%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure Settings %%%
%%%%%%%%%%%%%%%%%%%%%%%

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3)/1.6 scrsz(4)/1.18])%full screen: [1 scrsz(4)/2 scrsz(3)/1 scrsz(4)])
set(gcf,'paperpositionmode','auto'); 
TitleFont=16;
nbV=3;
nbH=length(TempTab);
spaceH=0.015;
spaceV=0.015;
borderH=0.15;
borderV=0.15;
HGraph = (1-(length(TempTab)-1)*spaceH-2*borderH)/3;
VGraph = (1-2*spaceV-2*borderV)/3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Normalisation and bounds setting %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%NormalisRate=1/(0.001*BinSize);
NormalisRate=1;
MaxFeature=max(max(FeaturesTab(:,BinIndex,TempTab,state)));
%colmap=jet(MaxFeature);
col = hot(MaxFeature+4);
colmap=col(1:MaxFeature,:);

for iy=1:length(TempTab)    
    nozero = find(PexpTab(:,BinIndex,TempTab(iy),state)>0);    
    xmin(iy) = min(min(PexpTab(nozero,BinIndex,TempTab(iy),state)));
    xmax(iy) = max(max(PexpTab(nozero,BinIndex,TempTab(iy),state)));

    ymin(iy) = min([min(PthTab(nozero,BinIndex,TempTab(iy),state)) min(PthNoJ1Tab(nozero,BinIndex,TempTab(iy),state)) min(PthNoJTab(nozero,BinIndex,TempTab(iy),state)) ]);
    ymax(iy) = max([max(max(PthTab(:,BinIndex,TempTab(iy),state))) max(max(PthNoJ1Tab(:,BinIndex,TempTab,state))) max(max(PthNoJTab(:,BinIndex,TempTab,state))) ]);

    xmin(iy) = 10^((log10(xmin(iy)))-0.2);
    ymin(iy) = 10^((log10(ymin(iy))));

    xmax(iy) = 10^((log10(xmax(iy)))+0.2);
    ymax(iy) = 10^((log10(ymax(iy)))+0.2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loop for plotting the matrix figure %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ix=1:3
    for iy=1:length(TempTab)
        
        % Loading the probability distribution
        Pexp = PexpTab(:,BinIndex,TempTab(iy),state);
        nozero = find(Pexp>0);
        Pexp = Pexp(nozero);
        Feature = FeaturesTab(:,BinIndex,TempTab(iy),state);
        Feature = Feature(nozero);
        if ix==1
            Pth = PthTab(:,BinIndex,TempTab(iy),state);
        end
        if ix==2
            Pth = PthNoJ1Tab(:,BinIndex,TempTab(iy),state);
        end
        if ix==3
            Pth = PthNoJTab(:,BinIndex,TempTab(iy),state);
        end        
        Pth = Pth(nozero);
        
        % Plotting the distribtion comparison 
        axes('position',[(borderV+(ix-1)*(VGraph+spaceV)) (borderH+((length(TempTab)+1-iy)-1)*(HGraph+spaceH)) VGraph HGraph]);
        for i=1:MaxFeature
            indF=find( Feature==i );
            if i==2
                hold on
            end
            loglog(NormalisRate*Pexp(indF),NormalisRate*Pth(indF),'.','Color',colmap(i,:),'MarkerSize',15)
            hold on
        end
        a=[min(NormalisRate*Pexp) max(NormalisRate*Pexp)];
        b=logspace(log10(min(Pexp)),log10(max(Pexp)),500);
        y1 = NormalisRate*b + NormalisRate*1.96.*sqrt(b.*(1-b)./(N_r-1));
        y2 = NormalisRate*b - NormalisRate*1.96.*sqrt(b.*(1-b)./(N_r-1));
        loglog(NormalisRate*b,y1,'k--','Linewidth',2)
        loglog(NormalisRate*b,y2,'k--','Linewidth',2)
        loglog(a,a,'k','Linewidth',2);        
        set(gca,'FontSize',10);
        set(gca,'Ylim',NormalisRate*[ymin(iy) ymax(iy)]);
        set(gca,'Xlim',NormalisRate*[xmin(iy) xmax(iy)]);

        % Using the specific panel setting
        if (ix==1)&&(iy==length(TempTab))
            xlabel('                       Observed pattern probability','FontSize',16);
            set(gca,'fontsize',12,'xtick',[1e-5,1e-3,1e-1],'ytick',[1e-9,1e-5,1e-1])
        end        
        if (ix==1)&&(iy==3)
            ylabel('                            Predicted pattern probability','FontSize',16);
        end
        if (ix>1)
            set(gca,'Ytick',[])
        end        
        if (iy<length(TempTab))
            set(gca,'Xtick',[])
        end
        if (ix==1)&&(iy==1)
            set(gca,'fontsize',12,'ytick',[1e-5,1e-3,1e-1])
        end
        if (ix==1)&&(iy==2)
            set(gca,'fontsize',12,'ytick',[1e-7,1e-4,1e-1])
        end
        if (ix==2)&&(iy==3)
            set(gca,'fontsize',12,'xtick',[1e-5,1e-3,1e-1])
        end
        if (ix==3)&&(iy==3)
            set(gca,'fontsize',12,'xtick',[1e-5,1e-3,1e-1])
        end
        if (iy==1)
            if (ix==1)
                title('Markov','FontSize',TitleFont);
            end
            if (ix==2)
                title('Ising','FontSize',TitleFont);
            end
            if (ix==3)
                title('Independent','FontSize',TitleFont);
            end
        end
    end
end

axes('Position',[0 0 1 1],'Visible','off');
h=colorbar('location','south');
set(h,'location','manual','position',[0.623 0.04 0.225 0.028])
colormap(colmap);
caxis([0 MaxFeature])
set(h,'fontsize',12,'xtick',(0.5:1:MaxFeature-0.5),'xticklabel',{(1:MaxFeature)})
text(0.74,0.09,'Spike number','HorizontalAlignment','center','FontSize',14)