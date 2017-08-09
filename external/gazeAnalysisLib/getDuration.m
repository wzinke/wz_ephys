function [dur] = getDuration(DATA, timecol)
% Function [dur] = getDuration(DATA, timecol)
%  Calculates the duration of the data-clip in the timeformat of timecol.

if isempty(DATA{1})
    dur = 0;
else
    dur = getValue(DATA, rowCount(DATA), timecol) - getValue(DATA, 1, timecol);
end