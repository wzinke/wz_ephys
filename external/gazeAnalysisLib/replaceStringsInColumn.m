function [DATA] = replaceStringsInColumn(DATA, column, value_before, value_after)
%Function [DATA] = replaceStringsInColumn(DATA, column, value_before, value_after)
%
% Replaces all the value_before's with value_after in the column specified
% by column (number) in DATA.

datavec = DATA{column};

datavec(strcmp(datavec, value_before)) = {value_after};

DATA{column} = datavec;