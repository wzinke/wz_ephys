function WIfunc_plotraster()

global handles

handles.nclasses = length(unique(handles.class_id));
classlabels = unique(handles.class_id);
linewidth=2;
colors = [0 0 1; ...
    1 0 0; ...
    0 .8 0; ...
    0 0 0.1724; ...
    1 0.1034 0.7241];
ntrials = length(handles.spiketimes);

if ~isfield(handles,'denoise')
    axes(handles.loadspikes.raster_subplot)
    cla
    hold on
    for trial_i = 1:length(handles.spiketimes)
        coloraux = ...
            find(handles.class_id(trial_i)==classlabels);
        classcolor=...
            colors(mod(coloraux-1,size(colors,1))+1,:);
        
        spiketimesaux = handles.spiketimes{trial_i};
        for spk_i = 1:length(spiketimesaux)
            plot([spiketimesaux(spk_i) spiketimesaux(spk_i)],...
                [trial_i-0.4 trial_i+0.4],...
                'color',classcolor,...
                'linewidth',linewidth)
        end
    end
else
    aux = handles.loadspikes.raster_subplot;
    set(get(aux,'children'),'color',[1 1 1]*.9)
    
    axes(handles.loadspikes.raster_subplot)
    for trial_i = 1:ntrials
        spiketimesaux = handles.denoise.denoised_spiketimes{trial_i};
        coloraux = ...
            find(handles.class_id(trial_i)==classlabels);
        classcolor=...
            colors(mod(coloraux-1,size(colors,1))+1,:);
        for spk_i = 1:length(spiketimesaux)
            plot([spiketimesaux(spk_i) spiketimesaux(spk_i)],...
                [trial_i-0.4 trial_i+0.4],...
                'color',classcolor,...
                'linewidth',linewidth)
        end
    end
    
end

lastspike = max(cellfun(@max,handles.spiketimes));
set(handles.loadspikes.raster_subplot,'xlim',[0 lastspike])
set(handles.loadspikes.raster_subplot,'ylim',[0 ntrials+1])

end

