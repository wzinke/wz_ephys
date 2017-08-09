function WIfunc_denoise()

global handles

handles.denoise.denoised_spiketimes = ...
    cell(1,length(handles.spiketimes));
handles.denoise.THRSdenoismask = ...
    handles.denoise.denoismask>handles.denoise.denoise_thrs;

lastspike = max(cellfun(@max,handles.spiketimes));
reclength = ...
    ceil(lastspike/handles.binsize)*handles.binsize;
bincenters = handles.binsize/2:handles.binsize:reclength;
ntrials = length(handles.spiketimes);

for trial_i = 1:ntrials
    spiketimesaux = handles.spiketimes{trial_i};
    for spk_i = 1:length(spiketimesaux)
        
        [~,bini] = min(abs(spiketimesaux(spk_i)-bincenters));
        
        if handles.denoise.THRSdenoismask(trial_i,bini)

            handles.denoise.denoised_spiketimes{trial_i} = ...
               [handles.denoise.denoised_spiketimes{trial_i} ...
               handles.spiketimes{trial_i}(spk_i)];
        end

    end
end


