function WIfunc_decode_leaveoneout()

global handles

if isfield(handles,'WIgui')
    set(handles.WIgui.maintext,'String',...
        'Decoding spike counts (leave-one out).')
end

spkcount = sum(handles.matrices.actmatrix,2);
class_labels = unique(handles.class_id);
nclasses = length(class_labels);
ntrials = length(spkcount);

dec_output = nan(ntrials,1);
for trial_i = 1:ntrials
    
    sample = spkcount(trial_i);
    trainingtrials = setdiff(1:ntrials,trial_i);    
    training = spkcount(trainingtrials);
    training = training+rand(size(training))/999999999;
    class_id_training = handles.class_id(trainingtrials);
    
    dec_output(trial_i) = ...
        classify(sample,training,class_id_training,'diagLinear');

end
    
    handles.decode.SPKCNTconfusionmatrix = ...
        confusionmat(handles.class_id,dec_output)';
    
    handles.decode.SPKCNTperf = ...
        sum(handles.class_id(:)==dec_output(:))/length(dec_output);
    
    if isfield(handles,'WIgui')
        set(handles.WIgui.maintext,'String',...
            'Decoding time patterns (Wavelets).')
    end
    
    for trial_i = 1:ntrials
        
        handles.decode.testingtrial = trial_i;
        handles.decode.trainingtrials = setdiff(1:ntrials,trial_i);
        WIfunc_select_coefs()
        
        sample = ...
            handles.matrices.wvmatrix(trial_i,...
            handles.matrices.selected_wvcoefs);
        
        training = ...
            handles.matrices.wvmatrix(handles.decode.trainingtrials,...
            handles.matrices.selected_wvcoefs);
        training = training+rand(size(training))/999999999;
        
        class_id_training = ...
            handles.class_id(handles.decode.trainingtrials);
        
        dec_output(trial_i) = ...
            classify(sample,training,class_id_training,'diagLinear');
        
    end
    
    handles.decode.WVconfusionmatrix = ...
        confusionmat(handles.class_id,dec_output)';

    handles.decode.WVperf = ...
        sum(handles.class_id(:)==dec_output(:))/length(dec_output);

end

