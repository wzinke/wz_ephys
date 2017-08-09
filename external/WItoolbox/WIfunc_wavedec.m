function WIfunc_wavedec()

global handles

ntrials = size(handles.matrices.actmatrix,1);

if exist('wavedec.m','file')

    trial_i = 1;
    [handles.matrices.wvmatrix,handles.matrices.L] = ...
        wavedec(handles.matrices.actmatrix(trial_i,:),...
        handles.nscales,'haar');
    handles.matrices.wvmatrix = ...
        repmat(handles.matrices.wvmatrix,ntrials,1);
    for trial_i = 2:ntrials
        handles.matrices.wvmatrix(trial_i,:) = ...
            wavedec(handles.matrices.actmatrix(trial_i,:),...
            handles.nscales,'haar');
    end
    
else
    
    trial_i = 1;
    [handles.matrices.wvmatrix,handles.matrices.L] = ...
        WIfunc_fix_wavedec(handles.matrices.actmatrix(trial_i,:),...
        handles.nscales);
    handles.matrices.wvmatrix = ...
        repmat(handles.matrices.wvmatrix,ntrials,1);
    for trial_i = 2:ntrials
        handles.matrices.wvmatrix(trial_i,:) = ...
            WIfunc_fix_wavedec(handles.matrices.actmatrix(trial_i,:),...
            handles.nscales);
    end
    
end

handles.matrices.wvmatrix_binsize = handles.binsize;
handles.matrices.wvmatrix_nscales = handles.nscales;



