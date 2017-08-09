function axis_MI = WIfunc_MI(response,class_id)
%%

nstim = length(unique(class_id));

axis_MI=NaN(size(response,2),1);
nbins_dist = 10;

for axis_i = 1:size(response,2)
    
    if var(response(:,axis_i))>0
        [~,histcounts] = hist(response(:,axis_i),nbins_dist);
    else
        [~,histcounts] = hist(0*response(:,axis_i),nbins_dist);
    end
    
    probabilites = NaN(nstim,nbins_dist);
    for stim_i = 1:nstim
        [probabilites(stim_i,:),~] = ...
            hist(response(class_id==stim_i,axis_i),...
            histcounts);
    end
      
    probabilites = probabilites./sum(sum(probabilites));
       
    MI_aux = 0;
    for bin_i = 1:nbins_dist
        p_value = sum(probabilites(:,bin_i)); % response probability
        for stim_i = 1:nstim
            p_stim = sum(probabilites(stim_i,:)); % stimulus probability
            if probabilites(stim_i,bin_i)>0
                PrGs = probabilites(stim_i,bin_i); % response probability given stim
                MI_aux = MI_aux + PrGs*log2(PrGs/(p_value*p_stim)); % MI
            end
        end
    end
    
    axis_MI(axis_i) = MI_aux;
    
end


