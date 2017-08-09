function WIfunc_decode_5050()

global handles

if isfield(handles,'WIgui')
    set(handles.WIgui.maintext,'String',...
        'Decoding spike counts.')
end

spkcount = sum(handles.matrices.actmatrix,2);
class_labels = unique(handles.class_id);
nclasses = length(class_labels);

handles.decode.trainingtrials = [];
handles.decode.testingtrials = [];
for class_i = 1:nclasses
    
    classlabel = class_labels(class_i);
    classtrials = find(handles.class_id==classlabel);
    nclasstrials = length(classtrials);
    
    [~,randorder] = sort(rand(nclasstrials,1));
    class_trainingtrials = ...
        classtrials(randorder(1:ceil(nclasstrials/2)));
    class_testingtrials = ...
        classtrials(randorder(1+ceil(nclasstrials/2):end));

    handles.decode.trainingtrials = ...
        [handles.decode.trainingtrials; class_trainingtrials(:)];
    handles.decode.testingtrials = ...
        [handles.decode.testingtrials; class_testingtrials(:)];

end

    sample = spkcount(handles.decode.testingtrials);
    class_id_sample = ...
        handles.class_id(handles.decode.testingtrials);
    training = spkcount(handles.decode.trainingtrials);
    training = training+rand(size(training))/999999999;
    class_id_training = ...
        handles.class_id(handles.decode.trainingtrials);
    
    dec_output = ...
        classify(sample,training,class_id_training,'diagLinear');
    handles.decode.SPKCNTconfusionmatrix = ...
        confusionmat(class_id_sample,dec_output);
    
    dec_output = ...
        classify(sample,training,class_id_training,'diagLinear');
    handles.decode.SPKCNTconfusionmatrix = ...
        confusionmat(class_id_sample,dec_output)';
    
    handles.decode.SPKCNTperf = ...
        sum(class_id_sample(:)==dec_output(:))/length(dec_output);
    
    WIfunc_select_coefs()
    
    if isfield(handles,'WIgui')
        set(handles.WIgui.maintext,'String',...
            'Decoding time patterns (Wavelets).')
    end
    
    sample = ...
        handles.matrices.wvmatrix(handles.decode.testingtrials,...
        handles.matrices.selected_wvcoefs);
    training = ...
        handles.matrices.wvmatrix(handles.decode.trainingtrials,...
        handles.matrices.selected_wvcoefs);
    training = training+rand(size(training))/999999999;
    dec_output = ...
        classify(sample,training,class_id_training,'diagLinear');
    handles.decode.WVconfusionmatrix = ...
        confusionmat(class_id_sample,dec_output)';

    handles.decode.WVperf = ...
        sum(class_id_sample(:)==dec_output(:))/length(dec_output);

end

