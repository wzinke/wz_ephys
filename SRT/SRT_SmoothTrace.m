function smTrace = SRT_SmoothTrace(Trace)
% Apply a filter to smooth the trace data.
%
% DESCRIPTION
%    Apply smoothing on the rows of <Trace>. 
%  
% SYNTAX
%
%   smTrace = SRT_SmoothTrace(Trace, span)
%
%   Input:  <Trace>  -  vector or matrix of continous data
%           
%
% .........................................................................
% wolf zinke, wolfzinke@gmail.com
%
% $Created : 19-Jun-2015 by wolf zinke
%

% ____________________________________________________________________________ %
%% 

% if(~exist('span','var') || isempty(span))
%     span = 100; % number of data points 
% end

if(isvector(Trace))
    smTrace = runSM(Trace);
else   
    smTrace = nan(size(Trace));
    for(t=1:size(Trace, 1))
        smTrace(t,:) = runSM(Trace(t,:));
    end
end
    
% ____________________________________________________________________________ %
function sV = runSM(V)

    sV = nan(size(V));
    
    ds = find(isfinite(V),1,'first');
    de = find(isfinite(V),1,'last');
    
    smdt = V(ds:de);
    
%     cs = s / length(smdt);
    
%     sV(ds:de) = smooth(smdt, cs,'loess');

    sV(ds:de) = smooth(smdt, 51,'sgolay', 16);
