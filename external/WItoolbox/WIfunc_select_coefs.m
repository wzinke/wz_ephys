function WIfunc_select_coefs()

global handles

MI_coefs = ...
    WIfunc_MI(handles.matrices.wvmatrix(handles.decode.trainingtrials,:),...
    handles.class_id(handles.decode.trainingtrials));
MI_thresholds = WIfunc_thresholdMI();

unbiased_MI = MI_coefs(:) - MI_thresholds(:);
handles.matrices.selected_wvcoefs = find(unbiased_MI>0);

if length(handles.matrices.selected_wvcoefs)>handles.maxwvcoefs
    [~,MIorder] = sort(unbiased_MI,'descend');
    handles.matrices.selected_wvcoefs = ...
        MIorder(1:handles.maxwvcoefs);
elseif length(handles.matrices.selected_wvcoefs)<handles.minwvcoefs
    [~,MIorder] = sort(unbiased_MI,'descend');
    handles.matrices.selected_wvcoefs = ...
        MIorder(1:handles.minwvcoefs);
end

if isfield(handles,'WIgui')
    set(handles.WIgui.maintext,'String',...
        'Wavelet coefficients selected.')
    
    
    if strcmp(handles.crossvalidation_method,'50/50')
        handles.decode.WVcoef_info_subplot = ...
            subplot('Position',[.065 .1 .3 .3]);
        cla
        nonselectedwvcs = setdiff(1:length(MI_coefs),...
            handles.matrices.selected_wvcoefs);
        plot(nonselectedwvcs,MI_coefs(nonselectedwvcs),'ok','linewidth',1), hold on
        plot(MI_thresholds,'--r','linewidth',3)
        plot(handles.matrices.selected_wvcoefs,...
            MI_coefs(handles.matrices.selected_wvcoefs),'or','linewidth',2),
        axis tight, box off
        set(handles.decode.WVcoef_info_subplot,...
            'fontsize',14)
        xlabel('Wavelet coefficient','fontsize',15)
        ylabel('Information (bits)','fontsize',15)
        set(handles.decode.WVcoef_info_subplot,'visible','on')
        pause(0.00001)
    end
end