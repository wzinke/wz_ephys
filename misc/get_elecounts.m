function [cnts, lbls] = get_elecounts(vec)
    cnts = [];
    lbls = sort(unique(vec));

    if(isnumeric(lbls))
        for(i=1:length(lbls))
            cnts(i) = sum(vec == lbls(i));
        end
    elseif(iscell(lbls))
        for(i=1:length(lbls))
            cnts(i) = sum(strcmp(vec,lbls(i)));
        end
    else
        error('counting elements of this data type not implemented yet!');
    end


