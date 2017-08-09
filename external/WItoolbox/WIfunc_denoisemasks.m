function WIfunc_denoisemasks

global handles

classlist = unique(handles.class_id);

handles.denoise.rec_meanPSTH = ...
    nan(handles.nclasses,size(handles.matrices.actmatrix,2));
for class_i = 1:handles.nclasses
    
    auxclass = handles.class_id == classlist(class_i);
    meanPSTH = ...
        sum(handles.matrices.actmatrix(auxclass,:));
    
    meanPSTH_wv = wavedec(meanPSTH,...
        handles.nscales,'haar');
    
    denoisd_meanPSTH_wv = ...
        zeros(size(meanPSTH_wv));
    denoisd_meanPSTH_wv(handles.matrices.selected_wvcoefs) = ...
        meanPSTH_wv(handles.matrices.selected_wvcoefs);
    
    handles.denoise.rec_meanPSTH(class_i,:) = ...
        waverec(denoisd_meanPSTH_wv,handles.matrices.L,'haar');
    
end

handles.denoise.denoise_thrs = ...
    handles.denoisthrs*std(abs(handles.denoise.rec_meanPSTH(:)));

handles.denoise.denoismask = ...
    nan(size(handles.matrices.actmatrix));
for class_i = 1:handles.nclasses
    
    auxclass = handles.class_id == classlist(class_i);
    auxmask = ...
        repmat(handles.denoise.rec_meanPSTH(class_i,:),...
        sum(auxclass),1);
    handles.denoise.denoismask(auxclass,:) = auxmask;
    
end