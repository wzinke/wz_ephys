function [meanrt, varrt, stert, fano, trialrate] = spk_get_meanresp(spktm, timwin, clip)

    windur = repmat(diff(timwin)+1,size(spktm,1),1);
    
    if(exist('clip','var') && ~isempty(clip))
        spktm(bsxfun(@gt, spktm, clip(:))) = NaN;
        
        windur2 = (clip(:) - timwin(1)) + 1;
        p = windur2 < windur;
        windur(p) = windur2(p);
    end
    
    pos = spktm >= timwin(1)  & spktm <= timwin(2);
    spkcnt    = sum(pos,2);
    trialrate = (1000 .* spkcnt) ./ windur;
    

    meanrt = nanmean(trialrate);
    varrt  =  nanvar(trialrate);
    stdrt  =  nanstd(trialrate);
    stert  =  stdrt / sqrt(sum(isfinite(trialrate)));
    fano   = varrt/meanrt;

