function ctrial = SPK_plot_spikeraster(spks, evtm, twin, startY, col, evcol)
    
    if(~exist('startY','var') || isempty(startY))
        startY = 0;
    end
    
    if(~exist('evtm','var'))
        evtm = [];
    end
    
    if(~exist('twin','var') || isempty(twin))
        twin = [min(spks(:)), max(spks(:))];
    end
    
    if(~exist('col','var') || isempty(col))
        col = [0, 0, 0];
    end
    
    if(~exist('evcol','var') || isempty(evcol))
        evcol = [1, 0, 0];
    end
            
    hold on;
    
    nTrials = size(spks,1);
    ctrial = startY;
    
    if(~isempty(evtm))
        med_evtm = nanmedian(evtm);
        plot([med_evtm,med_evtm], [startY,startY+nTrials+1], '-', 'color', evcol, 'LineWidth', 1.5);
        
        [~, Torder] = sort(evtm);
    else
        Torder = 1:nTrials;
    end
    
    if(~isempty(spks) && any(isfinite(spks(:))))

        for(s=1:nTrials)
            ctrial  = ctrial + 1;

            cspikes = spks(Torder(s),:);
            cspikes(isnan(cspikes)) = [];
            cspikes(cspikes<twin(1) | cspikes>twin(2)) = [];
            
            if(~isempty(evtm))
                plot([evtm(Torder(s)),evtm(Torder(s))], [ctrial-0.5,ctrial+0.5],'color',evcol, 'LineWidth', 1.5);
            end
            
            if(~isempty(cspikes))
                plot(cspikes,ctrial,'.', 'color', col);
            end

        end
    end

   