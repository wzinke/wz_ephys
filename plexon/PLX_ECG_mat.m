function [beatpeak] = PLX_ECG_mat(LFP, peakpos)

if(~exist('peakpos','var'))
    peakpos = [];
end

if(ischar(LFP))
    LFP = load(LFP) ;
end

beatpeak = {};

num_rows = 4;
num_cols = 6;
matsz    = 24;

figsz  = [0 0 num_cols*200 num_rows*150];
Ystart = linspace(97,3,num_rows+1)./100;
Ystart(1)   = [];
Xstart = linspace(3,97,num_cols+1)./100;
Xstart(end) = [];

xwd = min(diff(Xstart))-0.018;
ywd = min(abs(diff(Ystart)))-0.024;
figure('Position', figsz, 'Renderer', 'Painters');

chancnt = 0; % currently this assumes that channel count starts at 1. Needs to be corrected to be more flexible!
beatpeak = {};

for(x=1:num_cols)
    cX = Xstart(x);
    for(y=1:num_rows)
        cY = Ystart(y);
        
        chancnt = chancnt+1;
        
        if(chancnt > matsz)
            break
        end
        
        subplot('Position',[cX cY xwd ywd]);
        hold on;
        
        %         if(x>1)
        %             set(gca,'YTickLabel',[]);
        %         end
        
        if(y<num_rows)
            set(gca,'XTickLabel',[]);
        end
        
        beatpeak(chancnt) = {PLX_ECGavrg(LFP(chancnt,:), sprintf('LFP%02d%c', chancnt), peakpos)};
    end
end




