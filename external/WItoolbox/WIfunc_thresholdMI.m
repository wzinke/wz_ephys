function [ MI_thresholds ] = WIfunc_thresholdMI()

global handles

MI_coefs_surrogate = nan(handles.nsurr,size(handles.matrices.wvmatrix,2));
for surr_i = 1:handles.nsurr
    
    progressaux = num2str(round(100*(surr_i-1)/handles.nsurr));
    
    if isfield(handles,'WIgui')
        if strcmp(handles.crossvalidation_method,'50/50')
            set(handles.WIgui.maintext,'String',...
                ['Selecting wavelet coefficients. Progress: ' progressaux '%'])
            pause(.0000001)
        elseif strcmp(handles.crossvalidation_method,'Leave-one out')
            auxprint = ['Decoding trial ' num2str(handles.decode.testingtrial) '. '];
            set(handles.WIgui.maintext,'String',...
                [auxprint 'Selecting wavelet coefficients. Progress: ' progressaux '%'])
            pause(.0000001)
        end
    end
    
    class_id_surrogate = ...
        handles.class_id(randperm(length(handles.class_id)));
    MI_coefs_surrogate(surr_i,:) = ...
        WIfunc_MI(handles.matrices.wvmatrix,class_id_surrogate);
end

MI_thresholds = nan(1,size(handles.matrices.wvmatrix,2));
levelbounds = cumsum([0 handles.matrices.L(1:end-1)]);
levelthres = NaN(1,length(handles.matrices.L)-1);
for level_i=1:length(handles.matrices.L)-1
    level_coefs = levelbounds(level_i)+1:levelbounds(level_i+1);
    level_MI_nulldist = MI_coefs_surrogate(:,level_coefs);
    level_MI_nulldist = level_MI_nulldist(:);
    levelthres(level_i) = ...
        prctile(level_MI_nulldist,handles.percentile);
    MI_thresholds(level_coefs) = levelthres(level_i);
end


end

