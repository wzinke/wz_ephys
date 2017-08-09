function [value] = getValue(DATA, row, column)
%Function [value] = getValue(DATA, row, column)
%
% Returns a value specified by row and column from the array DATA.

raw_value = DATA{column}(row);

if iscell(raw_value)
    value = raw_value{1};
else
    value = raw_value;
end