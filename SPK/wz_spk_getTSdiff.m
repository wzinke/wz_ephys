function [TDTvec, TDTci] = wz_spk_getTSdiff(dvals, thrval, minseq, uBnd, lBnd, xdata)
% wz_spk_getTSdiff - identify the time point here a signal in a time series
%                    exceds a specified threshold for a sufficient number
%                    of subsequent time bins
%
% DESCRIPTION
%
% SYNTAX
%   TST = wz_spk_getTSdiff()
%
%   Input:
%
%
%   Output:
%
%
% REFERENCES
%
%
% .........................................................................
% wolf zinke, wolfzinke@gmail.com
%
% $Created : 06-Oct-2014 by wolf zinke
% $Modified:

TDTci = [NaN, NaN];

if(~exist('minseq','var') || isempty(minseq))
    minseq = 10;
end

if(~exist('xdata','var') || isempty(xdata))
    if(isvector(dvals))
        xdata = 1:length(dvals);
    else
        xdata = 1:size(dvals,2);
    end
end

% find bins exceeding the threshold
sigval  = dvals <= thrval;

for(r=1:size(sigval,1))
    % get the first bin that is followed by sufficient many super-threshold bins
    sigpos = strfind(sigval(r,:), ones(1,minseq));
    if(~isempty(sigpos))
        TDTvec(r) = xdata(sigpos(1));
    else
        TDTvec(r) = NaN;
    end
end

%% if dvals is a vector and upper and lower bounds are specified determine CI
if(length(TDTvec) ==1 && ~isnan(TDTvec(1)))
    
    if(exist('lBnd','var') && ~isempty(lBnd))
        if(isfinite(TDTest))
            sigpos = xdata(find(lBnd > thrval & xdata < TDTvec,1,'last')+1);
            if(~isempty(sigpos))
                TDTci(1) = xdata(sigpos);
            end
        end
    end
    
    if(exist('uBnd','var') && ~isempty(uBnd))
        if(isfinite(TDTest))
            sigpos = xdata(find(uBnd < thrval & xdata > TDTvec,1,'first'));
            if(~isempty(sigpos))
                TDTci(2) = xdata(sigpos);
            end
        end
    end
end

