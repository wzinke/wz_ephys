function Gaze = SRT_SaccadeDetection(EyeX, EyeY, StimOn, GoCue, xvec, span)
% Read in a plexon data file and organize it according to trials.
%
% DESCRIPTION
%
%
% SYNTAX
%
%   Gaze = SRT_SaccadeDetection(EyeX, EyeY)
%
%   Input:
%
% .........................................................................
% wolf zinke, wolfzinke@gmail.com
%
% $Created : 19-Jun-2015 by wolf zinke
%

% ____________________________________________________________________________ %
%% check input arguments

if(~exist('span','var') || isempty(span))
    span = 100; % number of data points
end

if(~exist('xvec','var') || isempty(xvec))
    if(isvector(EyeX))
        xvec = length(EyeX);
    else
        xvec = size(EyeX,2);
    end
end

if(any(size(EyeX) ~= size(EyeY)))
    error('X and Y traces have different dimensions!');
end

% ____________________________________________________________________________ %
%% smooth data (loess filter)

smX = SRT_SmoothTrace(EyeX, span);
smY = SRT_SmoothTrace(EyeY, span);

% get the distance relative to previous eye position
Dist = getDist(EyeX, EyeY);
smD  = getDist(smX, smY);

% ____________________________________________________________________________ %
%% determine relevant events
% get event position for stimulus onset
if(length(StimOn) > 1 && any(diff(StimOn)~=0))
    if(length(StimOn) ~=  size(EyeX, 1))
        error('Mismatch in number of StimOn times and number of trials!');
    end
    Son = bsxfun(@ge, xvec(:)', StimOn(:));
else
    Son = find(xvec >= StimOn(1), 1, 'first');
end

% get event position for go signal
if(length(GoCue) > 1 && any(diff(GoCue)~=0))
    if(length(GoCue) ~=  size(EyeX, 1))
        error('Mismatch in number of StimOn times and number of trials!');
    end
    Gon = bsxfun(@ge, xvec(:), GoCue(:));
else
    Gon = find(xvec >= GoCue(1), 1, 'first');
end


if(min(size(EyeX)) > 1) % trial matrix

else % single trace

end
% determine fixation median coordinate as reference

% get eccentricities relative to median fixation coordinate


% determine saccade start times (and points)

% determine saccade end times (and points)


% ____________________________________________________________________________ %
%% plot data
figure;




% ____________________________________________________________________________ %
%% helper functions
function D = getDist(X,Y)

    D = nan(size(X));
    D(:,2:end) = sqrt(diff(X,1,2).^2 + diff(Y,1,2).^2);

