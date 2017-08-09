function WIfunc_binmatrix() 

global handles

ntrials = length(handles.spiketimes);
lastspike = max(cellfun(@max,handles.spiketimes));
reclength = ...
    ceil(lastspike/handles.binsize)*handles.binsize;
bincenters = handles.binsize/2:handles.binsize:reclength;
nbins = length(bincenters);

handles.matrices.actmatrix = nan(ntrials,nbins);
for neuron_i = 1:ntrials
    handles.matrices.actmatrix(neuron_i,:) = ...
        hist(handles.spiketimes{neuron_i},bincenters);
end

handles.matrices.actmatrix_binsize = handles.binsize;




